#!/usr/bin/env python3
"""Run resumable NetBID2 QC for dynamically discovered PR67 sweep points."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any


DRIVERS = {
    "tf": ("BRCA100_TF.txt", "TF_", 2608),
    "sig": ("BRCA100_SIG.txt", "SIG_", 10680),
}
POINT_SCHEMA = "sjaracne-brca100-pr67-p-sweep-point-v1"
RUN_SCHEMA = "sjaracne-brca100-pr67-p-sweep-netbid2-v1"
SUMMARY_AGGREGATE_SCHEMA = (
    "sjaracne-brca100-pr67-p-sweep-netbid2-summary-aggregate-v1"
)
HTML_AGGREGATE_SCHEMA = (
    "sjaracne-brca100-pr67-p-sweep-netbid2-html-aggregate-v1"
)
REQUIRED_METRICS = {
    "candidate_drivers",
    "active_drivers",
    "active_driver_fraction",
    "edges",
    "incident_nodes",
    "weak_components",
    "largest_weak_component",
    "largest_weak_component_fraction",
    "density",
    "target_size_zero_mean",
    "target_size_zero_median",
    "target_size_zero_q25",
    "target_size_zero_q75",
    "target_size_zero_max",
    "target_size_active_mean",
    "target_size_active_median",
    "target_size_active_q25",
    "target_size_active_q75",
    "target_size_active_max",
    "scale_free_adjusted_r2",
}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprint(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    os.replace(partial, path)


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def discover_points(
    work_root: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any], str]:
    results_root = work_root / "results"
    sweep_design_path = work_root / "sweep_design.json"
    if not sweep_design_path.is_file():
        raise FileNotFoundError(sweep_design_path)
    sweep_design = load_json(sweep_design_path)
    if sweep_design.get("schema") != "sjaracne-brca100-pr67-p-sweep-v1":
        raise ValueError(f"Unexpected sweep design schema: {sweep_design_path}")
    design_points = sweep_design.get("all_points")
    if not isinstance(design_points, list) or not design_points:
        raise ValueError(f"Sweep design has no all_points list: {sweep_design_path}")
    if any(not isinstance(item, dict) or not isinstance(item.get("key"), str)
           for item in design_points):
        raise ValueError(f"Malformed all_points in {sweep_design_path}")
    design_by_key = {item["key"]: item for item in design_points}
    if len(design_by_key) != len(design_points):
        raise ValueError(f"Duplicate point key in {sweep_design_path}")

    manifest_paths = sorted(results_root.glob("*/point_manifest.json"))
    manifest_keys = {path.parent.name for path in manifest_paths}
    if manifest_keys != set(design_by_key):
        raise ValueError(
            "Point manifests do not exactly cover sweep_design.json all_points: "
            f"missing={sorted(set(design_by_key) - manifest_keys)}, "
            f"extra={sorted(manifest_keys - set(design_by_key))}"
        )
    points: list[dict[str, Any]] = []
    seen_probabilities: set[float] = set()
    for manifest_path in manifest_paths:
        manifest = load_json(manifest_path)
        key = manifest.get("key")
        if manifest.get("schema") != POINT_SCHEMA:
            raise ValueError(f"Unexpected point schema: {manifest_path}")
        if not isinstance(key, str) or key != manifest_path.parent.name:
            raise ValueError(f"Point key/path mismatch: {manifest_path}")
        design_point = design_by_key[key]
        mismatched_fields = [
            field for field, value in design_point.items()
            if manifest.get(field) != value
        ]
        if mismatched_fields:
            raise ValueError(
                f"Point manifest disagrees with sweep design at {manifest_path}: "
                + ", ".join(mismatched_fields)
            )
        probability = float(manifest["p_value"])
        cutoff = float(manifest["mi_cutoff"])
        if not math.isfinite(probability) or probability <= 0.0:
            raise ValueError(f"Invalid point probability: {manifest_path}")
        if not math.isfinite(cutoff) or cutoff <= 0.0:
            raise ValueError(f"Invalid point cutoff: {manifest_path}")
        if probability in seen_probabilities:
            raise ValueError(f"Duplicate point probability {probability:.17g}")
        seen_probabilities.add(probability)
        points.append(
            {
                "key": key,
                "p_value": probability,
                "mi_cutoff": cutoff,
                "manifest": manifest,
                "manifest_path": manifest_path,
                "root": manifest_path.parent,
            }
        )
    if not points:
        raise ValueError(f"No immutable point manifests found under {results_root}")
    return (
        sorted(points, key=lambda item: (item["p_value"], item["key"])),
        sweep_design,
        sha256_file(sweep_design_path),
    )


def select_values(
    specification: str, available: list[str], label: str, *, allow_none: bool = False
) -> list[str]:
    if specification == "all":
        return list(available)
    if allow_none and specification in ("", "none"):
        return []
    requested = [item.strip() for item in specification.split(",") if item.strip()]
    unknown = sorted(set(requested) - set(available))
    if unknown:
        raise ValueError(f"Unknown {label}(s): {', '.join(unknown)}")
    if not requested:
        raise ValueError(f"At least one {label} must be selected")
    if len(set(requested)) != len(requested):
        raise ValueError(f"Duplicate {label} selected")
    return requested


def probe_environment(wrapper: Path, r_script: Path) -> tuple[dict[str, str], str]:
    completed = subprocess.run(
        [str(wrapper), "Rscript", str(r_script), "--probe"],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "NetBID2 environment probe failed:\n" + completed.stderr.strip()
        )
    rows = list(csv.DictReader(completed.stdout.splitlines(), delimiter="\t"))
    if not rows or any(set(row) != {"component", "version"} for row in rows):
        raise ValueError("Unexpected NetBID2 environment probe output")
    environment = {row["component"]: row["version"] for row in rows}
    if (
        len(rows) != 4
        or len(environment) != len(rows)
        or set(environment) != {"R", "NetBID2", "NetBID2_remote_sha", "igraph"}
    ):
        raise ValueError(f"Incomplete NetBID2 environment probe: {environment}")
    return environment, completed.stderr


def inventory(root: Path) -> list[dict[str, object]]:
    if not root.is_dir():
        raise RuntimeError(f"Missing output directory: {root}")
    records = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.stat().st_size == 0:
            raise RuntimeError(f"Empty NetBID2 artifact: {path}")
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    if not records:
        raise RuntimeError(f"Empty NetBID2 output directory: {root}")
    return records


def read_network_summary(path: Path) -> dict[str, float]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or any(set(row) != {"metric", "value"} for row in rows):
        raise ValueError(f"Unexpected network summary format: {path}")
    if any(not row["metric"] for row in rows):
        raise ValueError(f"Empty metric in {path}")
    if len({row["metric"] for row in rows}) != len(rows):
        raise ValueError(f"Duplicate metric in {path}")
    metrics = {
        row["metric"]: float("nan") if row["value"] == "NA" else float(row["value"])
        for row in rows
    }
    if set(metrics) != REQUIRED_METRICS:
        raise ValueError(
            f"Unexpected metrics in {path}: missing={sorted(REQUIRED_METRICS - set(metrics))}, "
            f"extra={sorted(set(metrics) - REQUIRED_METRICS)}"
        )
    return metrics


def validate_output(
    root: Path,
    *,
    mode: str,
    prefix: str,
    driver_ids: list[str],
    expected_edges: int,
    expected_environment: dict[str, str],
) -> list[dict[str, object]]:
    required = [
        root / "network_summary.tsv",
        root / "driver_target_sizes.tsv",
        root / "netbid_environment.tsv",
    ]
    if mode == "html":
        required.append(root / f"{prefix}netQC.html")
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty NetBID2 output: {path}")

    metrics = read_network_summary(root / "network_summary.tsv")
    if int(metrics["candidate_drivers"]) != len(driver_ids):
        raise ValueError(f"Candidate-driver count mismatch in {root}")
    if int(metrics["edges"]) != expected_edges:
        raise ValueError(f"Consensus/NetBID2 edge-count mismatch in {root}")

    with (root / "driver_target_sizes.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if any(set(row) != {"driver", "target_count"} for row in rows):
        raise ValueError(f"Unexpected driver-target format in {root}")
    if [row["driver"] for row in rows] != driver_ids:
        raise ValueError(f"Driver order/content mismatch in {root}")
    target_counts = [int(row["target_count"]) for row in rows]
    if any(value < 0 for value in target_counts):
        raise ValueError(f"Negative driver target count in {root}")
    if sum(target_counts) != expected_edges:
        raise ValueError(f"Driver target counts do not sum to edge count in {root}")
    if sum(value > 0 for value in target_counts) != int(metrics["active_drivers"]):
        raise ValueError(f"Active-driver count mismatch in {root}")

    with (root / "netbid_environment.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        environment_rows = list(csv.DictReader(handle, delimiter="\t"))
    output_environment = {
        row.get("component"): row.get("version") for row in environment_rows
    }
    if (
        len(environment_rows) != len(expected_environment)
        or len(output_environment) != len(environment_rows)
        or any(set(row) != {"component", "version"} for row in environment_rows)
        or output_environment != expected_environment
    ):
        raise ValueError(f"Unexpected environment table in {root}")
    return inventory(root)


def validate_record(
    *,
    root: Path,
    record: dict[str, Any],
    expected_fingerprint: str,
    mode: str,
    prefix: str,
    driver_ids: list[str],
    expected_edges: int,
    stdout_path: Path,
    stderr_path: Path,
    expected_environment: dict[str, str],
) -> None:
    actual_inventory = validate_output(
        root,
        mode=mode,
        prefix=prefix,
        driver_ids=driver_ids,
        expected_edges=expected_edges,
        expected_environment=expected_environment,
    )
    for log_path in (stdout_path, stderr_path):
        if not log_path.is_file():
            raise RuntimeError(f"Missing NetBID2 log: {log_path}")
    expected = {
        "fingerprint": expected_fingerprint,
        "output_inventory": actual_inventory,
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
    }
    mismatches = [key for key, value in expected.items() if record.get(key) != value]
    if mismatches:
        raise RuntimeError(
            f"NetBID2 manifest validation failed at {root}: {', '.join(mismatches)}"
        )


def run_mode(
    *,
    point: dict[str, Any],
    driver: str,
    driver_file: Path,
    driver_ids: list[str],
    prefix: str,
    mode: str,
    wrapper: Path,
    wrapper_hash: str,
    r_script: Path,
    r_script_hash: str,
    environment: dict[str, str],
    sweep_design_hash: str,
) -> dict[str, Any]:
    arm_root = point["root"] / driver
    consensus = arm_root / "consensus/consensus_network_ncol_.txt"
    consensus_manifest_path = arm_root / "consensus_manifest.json"
    if not consensus.is_file() or not consensus_manifest_path.is_file():
        raise RuntimeError(f"Consensus is incomplete for {point['key']}/{driver}")
    consensus_manifest = load_json(consensus_manifest_path)
    consensus_hash = sha256_file(consensus)
    if (
        consensus_manifest.get("stage") != point["key"]
        or consensus_manifest.get("driver") != driver
        or consensus_manifest.get("ncol", {}).get("sha256") != consensus_hash
    ):
        raise RuntimeError(f"Consensus provenance mismatch for {point['key']}/{driver}")
    expected_edges = int(consensus_manifest["ncol"]["edges"])

    mode_suffix = "netbid2_qc" if mode == "summary" else "netbid2_qc_html"
    output_root = arm_root / mode_suffix
    partial_root = arm_root / f"{mode_suffix}.partial"
    manifest_path = arm_root / f"{mode_suffix}_manifest.json"
    pending_path = arm_root / f"{mode_suffix}_manifest.pending.json"
    log_root = arm_root / "logs"
    log_root.mkdir(parents=True, exist_ok=True)
    stdout_path = log_root / f"{mode_suffix}.stdout.log"
    stderr_path = log_root / f"{mode_suffix}.stderr.log"

    fingerprint_payload = {
        "schema": RUN_SCHEMA,
        "mode": mode,
        "point": point["key"],
        "p_value": point["p_value"],
        "mi_cutoff": point["mi_cutoff"],
        "point_manifest_sha256": sha256_file(point["manifest_path"]),
        "sweep_design_sha256": sweep_design_hash,
        "driver": driver,
        "driver_sha256": sha256_file(driver_file),
        "consensus_sha256": consensus_hash,
        "consensus_manifest_sha256": sha256_file(consensus_manifest_path),
        "r_script_sha256": r_script_hash,
        "wrapper_sha256": wrapper_hash,
        "environment": environment,
        "prefix": prefix,
    }
    run_fingerprint = fingerprint(fingerprint_payload)

    if manifest_path.is_file() and output_root.is_dir():
        record = load_json(manifest_path)
        validate_record(
            root=output_root,
            record=record,
            expected_fingerprint=run_fingerprint,
            mode=mode,
            prefix=prefix,
            driver_ids=driver_ids,
            expected_edges=expected_edges,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            expected_environment=environment,
        )
        if pending_path.exists():
            pending_record = load_json(pending_path)
            if pending_record != record:
                raise RuntimeError(
                    f"Pending/completed NetBID2 manifests disagree: {pending_path}"
                )
            pending_path.unlink()
        print(f"[NETBID2 {mode.upper()}] {point['key']}/{driver} resume", flush=True)
        return record

    if manifest_path.exists() != output_root.exists():
        if manifest_path.exists():
            raise RuntimeError(f"Manifest exists without output: {manifest_path}")
        if not pending_path.is_file():
            raise RuntimeError(f"Unverifiable orphan NetBID2 output: {output_root}")

    if pending_path.is_file() and (output_root.is_dir() or partial_root.is_dir()):
        if output_root.is_dir() and partial_root.exists():
            raise RuntimeError(f"Both final and partial NetBID2 outputs exist: {arm_root}")
        recovery_root = output_root if output_root.is_dir() else partial_root
        record = load_json(pending_path)
        validate_record(
            root=recovery_root,
            record=record,
            expected_fingerprint=run_fingerprint,
            mode=mode,
            prefix=prefix,
            driver_ids=driver_ids,
            expected_edges=expected_edges,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            expected_environment=environment,
        )
        if recovery_root == partial_root:
            os.replace(partial_root, output_root)
        record["recovered_after_interrupted_manifest"] = True
        # Keep pending and completed records byte-equivalent across the final
        # narrow crash window so a later M+O+P resume can verify before unlink.
        atomic_json(pending_path, record)
        atomic_json(manifest_path, record)
        pending_path.unlink()
        print(f"[NETBID2 {mode.upper()}] {point['key']}/{driver} recovered", flush=True)
        return record

    if output_root.exists() or manifest_path.exists() or pending_path.exists():
        raise RuntimeError(f"Unverifiable NetBID2 state under {arm_root}")
    if partial_root.exists():
        shutil.rmtree(partial_root)

    command = [
        str(wrapper),
        "Rscript",
        str(r_script),
        str(consensus),
        str(driver_file),
        str(partial_root),
        prefix,
        "true" if mode == "html" else "false",
    ]
    print(f"[NETBID2 {mode.upper()}] {point['key']}/{driver}", flush=True)
    with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, \
            stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
        subprocess.run(command, stdout=stdout, stderr=stderr, check=True)
    output_inventory = validate_output(
        partial_root,
        mode=mode,
        prefix=prefix,
        driver_ids=driver_ids,
        expected_edges=expected_edges,
        expected_environment=environment,
    )
    record = {
        **fingerprint_payload,
        "fingerprint": run_fingerprint,
        "command": command,
        "finished_at_utc": utc_now(),
        "output": str(output_root),
        "output_inventory": output_inventory,
        "stdout_sha256": sha256_file(stdout_path),
        "stderr_sha256": sha256_file(stderr_path),
        "stderr_bytes": stderr_path.stat().st_size,
    }
    atomic_json(pending_path, record)
    os.replace(partial_root, output_root)
    atomic_json(manifest_path, record)
    pending_path.unlink()
    return record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--benchmark-repo",
        type=Path,
        default=Path("/mnt/d/GitHub/SJARACNe-brca100-netbid-qc"),
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=(
            Path.home()
            / "sjaracne-benchmarks"
            / "brca100-pr67-threshold-sweep-20260818"
        ),
    )
    parser.add_argument("--points", default="all")
    parser.add_argument("--drivers", default="all")
    parser.add_argument(
        "--html-points",
        default="none",
        help="Comma-separated selected point keys, all, or none (default)",
    )
    return parser.parse_args()


def finalized_aggregate(payload: dict[str, Any]) -> dict[str, Any]:
    result = dict(payload)
    result["fingerprint"] = fingerprint(result)
    return result


def main() -> int:
    args = parse_args()
    results_root = args.work_root / "results"
    discovered, _sweep_design, sweep_design_hash = discover_points(args.work_root)
    point_map = {point["key"]: point for point in discovered}
    available_points = [point["key"] for point in discovered]
    point_keys = select_values(args.points, available_points, "point")
    driver_keys = select_values(args.drivers, list(DRIVERS), "driver")
    html_points = select_values(
        args.html_points, available_points, "HTML point", allow_none=True
    )
    outside_selection = sorted(set(html_points) - set(point_keys))
    if outside_selection:
        raise ValueError(
            "--html-points must be a subset of --points: "
            + ", ".join(outside_selection)
        )
    complete_selection = (
        point_keys == available_points and driver_keys == list(DRIVERS)
    )
    if html_points and not complete_selection:
        raise ValueError(
            "HTML generation requires --points all --drivers all so every "
            "stable summary record is revalidated first"
        )

    wrapper = args.benchmark_repo / "benchmarks/brca100_netbid_qc/netbid2-r"
    r_script = (
        args.benchmark_repo
        / "benchmarks/brca100_pr67_threshold_sweep/run_netbid_qc.R"
    )
    for required in (wrapper, r_script):
        if not required.is_file():
            raise FileNotFoundError(required)
    environment, probe_stderr = probe_environment(wrapper, r_script)
    if probe_stderr:
        print("[NETBID2 PROBE STDERR] " + probe_stderr.rstrip(), flush=True)
    wrapper_hash = sha256_file(wrapper)
    r_script_hash = sha256_file(r_script)

    summary_records: list[dict[str, Any]] = []
    summary_by_arm: dict[tuple[str, str], dict[str, Any]] = {}
    arm_context: dict[tuple[str, str], tuple[Path, list[str], str]] = {}
    for key in point_keys:
        point = point_map[key]
        inputs = point["manifest"].get("inputs", {})
        for driver in driver_keys:
            filename, prefix, expected_count = DRIVERS[driver]
            driver_file = args.work_root / "inputs" / filename
            if not driver_file.is_file():
                raise FileNotFoundError(driver_file)
            driver_ids = [
                line.strip()
                for line in driver_file.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            if len(driver_ids) != expected_count or len(set(driver_ids)) != len(driver_ids):
                raise ValueError(f"Unexpected {driver} driver list")
            expected_input_hash = inputs.get(filename, {}).get("sha256")
            if expected_input_hash != sha256_file(driver_file):
                raise RuntimeError(
                    f"Driver input does not match point manifest: {key}/{driver}"
                )
            summary = run_mode(
                point=point,
                driver=driver,
                driver_file=driver_file,
                driver_ids=driver_ids,
                prefix=prefix,
                mode="summary",
                wrapper=wrapper,
                wrapper_hash=wrapper_hash,
                r_script=r_script,
                r_script_hash=r_script_hash,
                environment=environment,
                sweep_design_hash=sweep_design_hash,
            )
            summary_records.append(summary)
            summary_by_arm[(key, driver)] = summary
            arm_context[(key, driver)] = (driver_file, driver_ids, prefix)

    summary_aggregate = finalized_aggregate({
        "schema": SUMMARY_AGGREGATE_SCHEMA,
        "environment": environment,
        "sweep_design_sha256": sweep_design_hash,
        "all_sweep_points": [item["key"] for item in discovered],
        "selection": {
            "points": point_keys,
            "drivers": driver_keys,
        },
        "summary_runs": summary_records,
    })
    summary_aggregate_path = results_root / "netbid2_qc_manifest.json"
    if html_points:
        if not summary_aggregate_path.is_file():
            raise RuntimeError(
                "Full summary aggregate is missing; first run with "
                "--points all --drivers all --html-points none"
            )
        existing_summary_aggregate = load_json(summary_aggregate_path)
        if existing_summary_aggregate != summary_aggregate:
            raise RuntimeError(
                "Stable NetBID2 summary aggregate is stale or incomplete; "
                "rerun summary-only QC before HTML generation"
            )
    elif complete_selection:
        if (
            not summary_aggregate_path.is_file()
            or load_json(summary_aggregate_path) != summary_aggregate
        ):
            atomic_json(summary_aggregate_path, summary_aggregate)
    else:
        print(
            "[NETBID2] partial summary selection: stable root aggregate unchanged",
            flush=True,
        )

    html_records: list[dict[str, Any]] = []
    for key in html_points:
        point = point_map[key]
        for driver in driver_keys:
            driver_file, driver_ids, prefix = arm_context[(key, driver)]
            html = run_mode(
                point=point,
                driver=driver,
                driver_file=driver_file,
                driver_ids=driver_ids,
                prefix=prefix,
                mode="html",
                wrapper=wrapper,
                wrapper_hash=wrapper_hash,
                r_script=r_script,
                r_script_hash=r_script_hash,
                environment=environment,
                sweep_design_hash=sweep_design_hash,
            )
            summary = summary_by_arm[(key, driver)]
            for filename_to_match in (
                "network_summary.tsv",
                "driver_target_sizes.tsv",
                "netbid_environment.tsv",
            ):
                summary_path = Path(summary["output"]) / filename_to_match
                html_path = Path(html["output"]) / filename_to_match
                if sha256_file(summary_path) != sha256_file(html_path):
                    raise RuntimeError(
                        f"Summary/full HTML mismatch: {key}/{driver}/"
                        f"{filename_to_match}"
                    )
            html_records.append(html)

    if html_points:
        html_aggregate = finalized_aggregate({
            "schema": HTML_AGGREGATE_SCHEMA,
            "environment": environment,
            "sweep_design_sha256": sweep_design_hash,
            "all_sweep_points": [item["key"] for item in discovered],
            "selection": {
                "points": point_keys,
                "drivers": driver_keys,
                "html_points": html_points,
            },
            "html_runs": html_records,
        })
        atomic_json(
            results_root / "netbid2_qc_html_manifest.json", html_aggregate
        )
    print(
        f"[NETBID2] complete: {len(summary_records)} summary arm(s), "
        f"{len(html_records)} HTML arm(s)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
