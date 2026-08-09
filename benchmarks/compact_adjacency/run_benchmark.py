#!/usr/bin/env python3

"""Benchmark compact sorted adjacency rows against the map-row baseline."""

import argparse
import csv
import hashlib
import json
import os
import statistics
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


TIME_MARKER = "__SJARACNE_GNU_TIME__"


@dataclass(frozen=True)
class BenchmarkCase:
    case_id: str
    suite: str
    mode: str
    expression: Path
    adjacency: Path | None
    hubs: Path | None
    tf_list: Path | None
    extra_args: tuple[str, ...]
    source_rows: str = ""
    genes_in_graph: str = ""
    directed_edges: str = ""


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run serial, matched map-row/compact-row benchmarks under GNU time."
        )
    )
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--synthetic-expression", required=True, type=Path)
    parser.add_argument("--fixtures-dir", required=True, type=Path)
    parser.add_argument("--brca-expression", required=True, type=Path)
    parser.add_argument("--brca-hubs", required=True, type=Path)
    parser.add_argument("--brca-adjacency", required=True, type=Path)
    parser.add_argument("--config-dir", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--warmups",
        type=int,
        default=1,
        help="Warmups per implementation for synthetic and BRCA replay cases.",
    )
    parser.add_argument(
        "--full-warmups",
        type=int,
        default=0,
        help="Warmups per implementation for the expensive full BRCA run.",
    )
    parser.add_argument("--time-command", type=Path, default=Path("/usr/bin/time"))
    parser.add_argument("--keep-outputs", action="store_true")
    return parser.parse_args()


def require_file(path, description):
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")
    return path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def network_data_sha256(path):
    """Hash non-header network records, preserving byte-level record order."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for line in handle:
            if line.strip() and not line.startswith(b">"):
                digest.update(line)
    return digest.hexdigest()


def load_cases(args):
    fixtures_dir = args.fixtures_dir.resolve()
    manifest = require_file(fixtures_dir / "cases.csv", "synthetic manifest")
    expression = require_file(args.synthetic_expression, "synthetic expression")

    cases = []
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            adjacency = require_file(
                fixtures_dir / row["adjacency"],
                f"adjacency for {row['case_id']}",
            )
            hubs = None
            if row["hubs"]:
                hubs = require_file(
                    fixtures_dir / row["hubs"], f"hubs for {row['case_id']}"
                )
            cases.append(
                BenchmarkCase(
                    case_id=row["case_id"],
                    suite="synthetic_replay",
                    mode=row["mode"],
                    expression=expression,
                    adjacency=adjacency,
                    hubs=hubs,
                    tf_list=None,
                    extra_args=("-p", "1", "-e", "0"),
                    source_rows=row["source_rows"],
                    genes_in_graph=row["genes_in_graph"],
                    directed_edges=row["directed_edges"],
                )
            )

    brca_expression = require_file(args.brca_expression, "BRCA100 expression")
    brca_hubs = require_file(args.brca_hubs, "BRCA100 hub/TF list")
    brca_adjacency = require_file(args.brca_adjacency, "BRCA100 adjacency")
    config_dir = args.config_dir.resolve()
    require_file(config_dir / "config_threshold.txt", "MI threshold configuration")
    config_arg = str(config_dir) + os.sep
    brca_args = (
        "-p",
        "1e-7",
        "-e",
        "0",
        "-a",
        "adaptive_partitioning",
        "-r",
        "1",
        "-H",
        config_arg,
        "-N",
        "40",
        "-S",
        "1",
    )
    cases.extend(
        (
            BenchmarkCase(
                case_id="brca100_adjacency_replay",
                suite="brca_replay",
                mode="selected_hubs_imported_adjacency",
                expression=brca_expression,
                adjacency=brca_adjacency,
                hubs=brca_hubs,
                tf_list=brca_hubs,
                # The imported file already contains retained MI edges. p=1
                # sets the load threshold to zero so this replay measures the
                # complete stored topology rather than filtering it again.
                extra_args=("-p", "1", "-e", "0"),
            ),
            BenchmarkCase(
                case_id="brca100_full_seed1",
                suite="brca_full",
                mode="selected_hubs_mi_and_dpi",
                expression=brca_expression,
                adjacency=None,
                hubs=brca_hubs,
                tf_list=brca_hubs,
                extra_args=brca_args,
            ),
        )
    )
    return cases


def build_command(binary, case, output):
    command = [str(binary), "-i", str(case.expression)]
    if case.adjacency is not None:
        command.extend(("-j", str(case.adjacency)))
    if case.hubs is not None:
        command.extend(("-s", str(case.hubs)))
    if case.tf_list is not None:
        command.extend(("-l", str(case.tf_list)))
    command.extend(case.extra_args)
    command.extend(("-o", str(output)))
    return command


def parse_time_metrics(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    marked = [line for line in lines if line.startswith(TIME_MARKER + "\t")]
    if len(marked) != 1:
        raise RuntimeError(f"Could not parse GNU time output {path}: {lines!r}")
    _, elapsed, user, system, max_rss = marked[0].split("\t")
    return {
        "gnu_elapsed_seconds": float(elapsed),
        "gnu_user_seconds": float(user),
        "gnu_system_seconds": float(system),
        "gnu_max_rss_kib": int(max_rss),
    }


def execute(
    *,
    binary,
    case,
    implementation,
    run_label,
    order_position,
    results_dir,
    time_command,
    keep_output,
):
    outputs_dir = results_dir / "outputs"
    logs_dir = results_dir / "logs"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    stem = f"{case.case_id}_{implementation}_{run_label}"
    output = outputs_dir / f"{stem}.adj"
    time_path = logs_dir / f"{stem}.time.tsv"
    stdout_path = logs_dir / f"{stem}.stdout.txt"
    stderr_path = logs_dir / f"{stem}.stderr.txt"
    command = build_command(binary, case, output)
    timed_command = [
        str(time_command),
        "-f",
        f"{TIME_MARKER}\t%e\t%U\t%S\t%M",
        "-o",
        str(time_path),
        "--",
        *command,
    ]

    started = time.perf_counter()
    result = subprocess.run(timed_command, capture_output=True, text=True)
    perf_wall = time.perf_counter() - started
    stdout_path.write_text(result.stdout, encoding="utf-8")
    stderr_path.write_text(result.stderr, encoding="utf-8")

    if result.returncode != 0:
        raise RuntimeError(
            f"{case.case_id} {implementation} {run_label} failed with "
            f"exit code {result.returncode}; see {stdout_path} and {stderr_path}"
        )
    if "[NETWORK] Applying DPI" not in result.stdout:
        raise RuntimeError(
            f"DPI did not run for {case.case_id} {implementation} {run_label}"
        )
    if not output.is_file():
        raise RuntimeError(f"Expected network output was not created: {output}")

    metrics = parse_time_metrics(time_path)
    record = {
        "case_id": case.case_id,
        "suite": case.suite,
        "mode": case.mode,
        "source_rows": case.source_rows,
        "genes_in_graph": case.genes_in_graph,
        "directed_edges": case.directed_edges,
        "implementation": implementation,
        "run_label": run_label,
        "order_position": order_position,
        "perf_wall_seconds": perf_wall,
        **metrics,
        "output_bytes": output.stat().st_size,
        "file_sha256": sha256(output),
        "network_data_sha256": network_data_sha256(output),
        "command": json.dumps(command),
    }
    if not keep_output:
        output.unlink()
    return record


def implementation_order(case_index, run_index):
    order = ("baseline", "candidate")
    if (case_index + run_index) % 2:
        order = tuple(reversed(order))
    return order


def write_csv(path, records):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def summarize(cases, records):
    summaries = []
    for case in cases:
        case_records = [row for row in records if row["case_id"] == case.case_id]
        hashes = {row["network_data_sha256"] for row in case_records}
        if len(hashes) != 1:
            raise RuntimeError(
                f"Network outputs differ for {case.case_id}: {sorted(hashes)}"
            )

        by_implementation = {
            name: [row for row in case_records if row["implementation"] == name]
            for name in ("baseline", "candidate")
        }

        def median(name, metric):
            return statistics.median(row[metric] for row in by_implementation[name])

        baseline_wall = median("baseline", "perf_wall_seconds")
        candidate_wall = median("candidate", "perf_wall_seconds")
        baseline_user = median("baseline", "gnu_user_seconds")
        candidate_user = median("candidate", "gnu_user_seconds")
        baseline_rss = median("baseline", "gnu_max_rss_kib")
        candidate_rss = median("candidate", "gnu_max_rss_kib")
        summaries.append(
            {
                "case_id": case.case_id,
                "suite": case.suite,
                "mode": case.mode,
                "source_rows": case.source_rows,
                "genes_in_graph": case.genes_in_graph,
                "directed_edges": case.directed_edges,
                "baseline_median_perf_wall_seconds": f"{baseline_wall:.6f}",
                "candidate_median_perf_wall_seconds": f"{candidate_wall:.6f}",
                "wall_speedup": f"{baseline_wall / candidate_wall:.6f}",
                "baseline_median_gnu_user_seconds": f"{baseline_user:.6f}",
                "candidate_median_gnu_user_seconds": f"{candidate_user:.6f}",
                "user_speedup": f"{baseline_user / candidate_user:.6f}",
                "baseline_median_gnu_system_seconds": (
                    f"{median('baseline', 'gnu_system_seconds'):.6f}"
                ),
                "candidate_median_gnu_system_seconds": (
                    f"{median('candidate', 'gnu_system_seconds'):.6f}"
                ),
                "baseline_median_max_rss_kib": f"{baseline_rss:.0f}",
                "candidate_median_max_rss_kib": f"{candidate_rss:.0f}",
                "candidate_rss_ratio": f"{candidate_rss / baseline_rss:.6f}",
                "rss_reduction_percent": (
                    f"{100.0 * (baseline_rss - candidate_rss) / baseline_rss:.6f}"
                ),
                "network_data_sha256": next(iter(hashes)),
            }
        )
    return summaries


def main():
    args = parse_args()
    if args.repeats < 1 or args.warmups < 0 or args.full_warmups < 0:
        raise ValueError("repeats must be positive and warmup counts non-negative")

    baseline = require_file(args.baseline, "baseline binary")
    candidate = require_file(args.candidate, "candidate binary")
    time_command = require_file(args.time_command, "GNU time executable")
    implementations = {"baseline": baseline, "candidate": candidate}
    cases = load_cases(args)
    results_dir = args.results_dir.resolve()
    results_dir.mkdir(parents=True, exist_ok=True)

    metadata = {
        "baseline": str(baseline),
        "baseline_sha256": sha256(baseline),
        "candidate": str(candidate),
        "candidate_sha256": sha256(candidate),
        "time_command": str(time_command),
        "repeats": args.repeats,
        "warmups": args.warmups,
        "full_warmups": args.full_warmups,
        "input_sha256": {
            "synthetic_expression": sha256(args.synthetic_expression.resolve()),
            "synthetic_manifest": sha256(
                (args.fixtures_dir.resolve() / "cases.csv")
            ),
            "brca_expression": sha256(args.brca_expression.resolve()),
            "brca_hubs": sha256(args.brca_hubs.resolve()),
            "brca_adjacency": sha256(args.brca_adjacency.resolve()),
            "config_threshold": sha256(
                args.config_dir.resolve() / "config_threshold.txt"
            ),
        },
    }
    (results_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )

    records = []
    for case_index, case in enumerate(cases):
        warmup_count = args.full_warmups if case.suite == "brca_full" else args.warmups
        for warmup in range(1, warmup_count + 1):
            for order_position, implementation in enumerate(
                implementation_order(case_index, warmup), start=1
            ):
                print(
                    f"warmup {case.case_id} {implementation} {warmup}/{warmup_count}",
                    flush=True,
                )
                execute(
                    binary=implementations[implementation],
                    case=case,
                    implementation=implementation,
                    run_label=f"warmup{warmup}",
                    order_position=order_position,
                    results_dir=results_dir,
                    time_command=time_command,
                    keep_output=args.keep_outputs,
                )

        for repeat in range(1, args.repeats + 1):
            for order_position, implementation in enumerate(
                implementation_order(case_index, repeat), start=1
            ):
                record = execute(
                    binary=implementations[implementation],
                    case=case,
                    implementation=implementation,
                    run_label=f"repeat{repeat}",
                    order_position=order_position,
                    results_dir=results_dir,
                    time_command=time_command,
                    keep_output=args.keep_outputs,
                )
                record["repeat"] = repeat
                records.append(record)
                print(
                    f"{case.case_id} {implementation} r{repeat}: "
                    f"wall={record['perf_wall_seconds']:.3f}s, "
                    f"max_rss={record['gnu_max_rss_kib']} KiB",
                    flush=True,
                )

    write_csv(results_dir / "raw_timings.csv", records)
    write_csv(results_dir / "summary.csv", summarize(cases, records))


if __name__ == "__main__":
    main()
