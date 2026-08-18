#!/usr/bin/env python3
"""Summarize a matched PR67 per-subsample p-value threshold sweep.

The analysis is deliberately limited to network/statistical QC.  It does not
select an operating threshold and does not claim biological validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "sjaracne-brca100-pr67-p-sweep-v1"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DRIVERS = ("tf", "sig")
CANDIDATE_FILES = {"tf": "BRCA100_TF.txt", "sig": "BRCA100_SIG.txt"}
NETWORK_FILENAME = "consensus_network_ncol_.txt"
SUPPORT_FILENAME = "consensus_support.tsv"
SVG_METADATA = {
    "Creator": "SJARACNe BRCA100 PR67 threshold sweep",
    "Date": None,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def first_present(record: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in record and record[key] is not None:
            return record[key]
    return None


def finite_float(value: Any, *, field: str, source: Path) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {field} in {source}: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {field} in {source}: {value!r}")
    return result


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read JSON object {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def point_metadata(results_root: Path) -> list[dict[str, Any]]:
    """Discover sweep points; point_manifest.json is the p-value authority."""
    points: list[dict[str, Any]] = []
    for point_root in sorted(path for path in results_root.iterdir() if path.is_dir()):
        manifest_path = point_root / "point_manifest.json"
        if not manifest_path.is_file():
            continue
        record = read_json(manifest_path)
        p_value = finite_float(record.get("p_value"), field="p_value", source=manifest_path)
        if not 0.0 < p_value < 1.0:
            raise ValueError(f"p_value must be between zero and one in {manifest_path}")
        for driver in DRIVERS:
            if not (point_root / driver).is_dir():
                raise ValueError(f"Missing {driver} arm under {point_root}")
        mi_cutoff_value = first_present(record, "mi_cutoff", "MI_cutoff", "cutoff")
        mi_cutoff = (
            math.nan
            if mi_cutoff_value is None
            else finite_float(mi_cutoff_value, field="mi_cutoff", source=manifest_path)
        )
        points.append(
            {
                "p_key": point_root.name,
                "p_value": p_value,
                "log10_p": math.log10(p_value),
                "p_token": str(record.get("p_token", point_root.name)),
                "p_label": str(
                    record.get("label", record.get("p_token", point_root.name))
                ),
                "role": str(record.get("role", "")),
                "mi_cutoff": mi_cutoff,
                "validation_class": str(record.get("validation_class", "")),
                "tail_extrapolated": bool(record.get("tail_extrapolated", False)),
                "commit": str(first_present(record, "commit", "commit_sha") or ""),
                "model_sha256": str(
                    first_present(record, "model_sha256", "null_model_sha256") or ""
                ),
                "manifest_sha256": sha256_file(manifest_path),
                "manifest": record,
                "root": point_root,
            }
        )
    if len(points) < 2:
        raise ValueError(
            f"Expected at least two results/<p_key>/point_manifest.json files under "
            f"{results_root}; found {len(points)}"
        )
    points.sort(key=lambda point: (point["p_value"], point["p_key"]))
    values = [point["p_value"] for point in points]
    if len(values) != len(set(values)):
        raise ValueError(f"Duplicate p_value entries in {results_root}")
    for invariant in ("commit", "model_sha256"):
        observed = {point[invariant] for point in points if point[invariant]}
        if len(observed) > 1:
            raise ValueError(
                f"The sweep changes {invariant} as well as p_value: {sorted(observed)}"
            )
    return points


def candidate_lists(input_root: Path) -> dict[str, list[str]]:
    candidates: dict[str, list[str]] = {}
    for driver, filename in CANDIDATE_FILES.items():
        path = input_root / filename
        if not path.is_file():
            raise ValueError(f"Missing fixed candidate list: {path}")
        values = [
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if not values or len(values) != len(set(values)):
            raise ValueError(f"Empty or duplicate-containing candidate list: {path}")
        candidates[driver] = values
    return candidates


def expression_identifiers(path: Path) -> set[str]:
    if not path.is_file():
        raise ValueError(f"Missing fixed expression matrix: {path}")
    identifiers: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        header = handle.readline()
        if not header or "\t" not in header:
            raise ValueError(f"Invalid expression-matrix header: {path}")
        for line_number, line in enumerate(handle, 2):
            identifier = line.split("\t", 1)[0].strip()
            if not identifier or identifier in identifiers:
                raise ValueError(
                    f"Empty or duplicate expression ID at {path}:{line_number}"
                )
            identifiers.add(identifier)
    if not identifiers:
        raise ValueError(f"Expression matrix contains no data rows: {path}")
    return identifiers


def arm_file(arm_root: Path, filename: str) -> Path:
    candidates = (arm_root / "consensus" / filename, arm_root / filename)
    matches = [path for path in candidates if path.is_file()]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {filename} directly under {arm_root} or its "
            f"consensus directory; found {matches}"
        )
    return matches[0]


def load_network(arm_root: Path) -> pd.DataFrame:
    path = arm_file(arm_root, NETWORK_FILENAME)
    frame = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
    required = {"source", "target", "MI"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing {sorted(required - set(frame.columns))} in {path}")
    if frame[["source", "target"]].isna().any().any():
        raise ValueError(f"Missing edge endpoint in {path}")
    if frame.duplicated(["source", "target"]).any():
        raise ValueError(f"Duplicate directed edges in {path}")
    frame["MI"] = pd.to_numeric(frame["MI"], errors="raise")
    if len(frame) and (not np.isfinite(frame["MI"]).all() or (frame["MI"] <= 0).any()):
        raise ValueError(f"MI must be finite and positive in {path}")
    return frame.set_index(["source", "target"], drop=False).sort_index()


def load_support(arm_root: Path) -> pd.DataFrame:
    path = arm_file(arm_root, SUPPORT_FILENAME)
    frame = pd.read_csv(path, sep="\t", dtype={"source": str, "target": str})
    required = {"source", "target", "support_count", "support_fraction"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing {sorted(required - set(frame.columns))} in {path}")
    if frame[["source", "target"]].isna().any().any():
        raise ValueError(f"Missing support edge endpoint in {path}")
    if frame.duplicated(["source", "target"]).any():
        raise ValueError(f"Duplicate directed edges in {path}")
    for column in ("support_count", "support_fraction"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if len(frame) and (
        not np.isfinite(frame[["support_count", "support_fraction"]]).all().all()
        or (frame["support_count"] < 1).any()
        or not frame["support_fraction"].between(0.0, 1.0).all()
    ):
        raise ValueError(f"Invalid support values in {path}")
    if "consensus_MI_roundtrip_match" in frame.columns and len(frame):
        if not (pd.to_numeric(frame["consensus_MI_roundtrip_match"]) == 1).all():
            raise ValueError(f"Consensus MI roundtrip mismatch in {path}")
    return frame.set_index(["source", "target"], drop=False).sort_index()


def merge_network_support(arm_root: Path) -> pd.DataFrame:
    network = load_network(arm_root)
    support = load_support(arm_root)
    if not network.index.equals(support.index):
        network_edges = set(network.index)
        support_edges = set(support.index)
        raise ValueError(
            f"Network/support edge mismatch in {arm_root}: "
            f"network-only={len(network_edges - support_edges)}, "
            f"support-only={len(support_edges - network_edges)}"
        )
    result = network.copy()
    result["support_count"] = support["support_count"]
    result["support_fraction"] = support["support_fraction"]
    if "consensus_MI" in support.columns and len(result):
        support_mi = pd.to_numeric(support["consensus_MI"], errors="raise")
        if not np.allclose(result["MI"], support_mi, rtol=0.0, atol=5.00001e-5):
            raise ValueError(f"Network/support MI mismatch in {arm_root}")
    return result


def command_p_value(record: dict[str, Any]) -> float | None:
    command = record.get("command")
    if not isinstance(command, list):
        return None
    try:
        position = command.index("-p")
    except ValueError:
        return None
    if position + 1 >= len(command):
        return None
    try:
        return float(command[position + 1])
    except (TypeError, ValueError):
        return None


def seeds_from_metadata(arm_root: Path, expected_p: float) -> tuple[set[int], str, bool]:
    metadata_root = arm_root / "seed_metadata"
    paths = sorted(metadata_root.glob("*.json")) if metadata_root.is_dir() else []
    if not paths:
        return set(), "", False
    seeds: set[int] = set()
    checked_commands = 0
    for path in paths:
        record = read_json(path)
        try:
            seed = int(record["seed"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"Missing or invalid seed in {path}") from error
        if seed in seeds:
            raise ValueError(f"Duplicate seed {seed} under {metadata_root}")
        seeds.add(seed)
        observed_p = command_p_value(record)
        if observed_p is not None:
            checked_commands += 1
            if not math.isclose(observed_p, expected_p, rel_tol=1e-12, abs_tol=0.0):
                raise ValueError(
                    f"Seed command p={observed_p:g} disagrees with point-manifest "
                    f"p={expected_p:g} in {path}"
                )
    return seeds, str(metadata_root), checked_commands == len(paths)


def aggregate_manifest_rows(
    manifest_path: Path,
    point_key: str,
    driver: str,
    expected_p: float,
) -> tuple[set[int], bool]:
    if not manifest_path.is_file():
        return set(), False
    frame = pd.read_csv(manifest_path, sep="\t", dtype=str)
    driver_column = next((name for name in ("driver", "driver_class") if name in frame), None)
    seed_column = "seed" if "seed" in frame else None
    key_column = next(
        (
            name
            for name in ("p_key", "point", "stage", "arm", "threshold_key")
            if name in frame
        ),
        None,
    )
    if driver_column is None or seed_column is None:
        return set(), False
    selected = frame[frame[driver_column] == driver]
    if key_column is not None:
        selected = selected[selected[key_column] == point_key]
    if selected.empty:
        return set(), False
    seeds = pd.to_numeric(selected[seed_column], errors="raise").astype(int)
    if seeds.duplicated().any():
        raise ValueError(f"Duplicate seeds for {point_key}/{driver} in {manifest_path}")
    p_column = next(
        (name for name in ("bootstrap_p", "p_value", "per_subsample_p") if name in selected),
        None,
    )
    p_checked = False
    if p_column is not None:
        observed = pd.to_numeric(selected[p_column], errors="raise").unique()
        if len(observed) != 1 or not math.isclose(
            float(observed[0]), expected_p, rel_tol=1e-12, abs_tol=0.0
        ):
            raise ValueError(
                f"Manifest p values disagree with {point_key} p={expected_p:g}: "
                f"{observed.tolist()} in {manifest_path}"
            )
        p_checked = True
    return set(seeds), p_checked


def collect_seed_sets(
    work_root: Path,
    points: list[dict[str, Any]],
) -> tuple[dict[tuple[str, str], set[int]], pd.DataFrame]:
    aggregate_paths = (
        work_root / "results" / "run_manifest.tsv",
        work_root / "run_manifest.tsv",
    )
    seed_sets: dict[tuple[str, str], set[int]] = {}
    summary_rows: list[dict[str, Any]] = []
    for point in points:
        for driver in DRIVERS:
            arm_root = point["root"] / driver
            seeds: set[int] = set()
            source = ""
            p_checked = False
            manifest_paths = (
                *aggregate_paths,
                point["root"] / "run_manifest.tsv",
                arm_root / "run_manifest.tsv",
            )
            for path in manifest_paths:
                manifest_seeds, manifest_p_checked = aggregate_manifest_rows(
                    path, point["p_key"], driver, point["p_value"]
                )
                if manifest_seeds:
                    seeds = manifest_seeds
                    source = str(path)
                    p_checked = manifest_p_checked
                    break
            if not seeds:
                seeds, source, p_checked = seeds_from_metadata(
                    arm_root, point["p_value"]
                )
            if not seeds:
                manifest_seeds = point["manifest"].get("seeds")
                if isinstance(manifest_seeds, list) and manifest_seeds:
                    seeds = {int(seed) for seed in manifest_seeds}
                    source = str(point["root"] / "point_manifest.json") + ":seeds"
            if not seeds:
                raise ValueError(
                    f"No seed evidence for {point['p_key']}/{driver}; expected an "
                    "aggregated/per-point run_manifest.tsv, seed_metadata/*.json, "
                    "or point_manifest.json seeds"
                )
            key = (point["p_key"], driver)
            seed_sets[key] = seeds
            encoded = ",".join(str(seed) for seed in sorted(seeds)).encode("ascii")
            summary_rows.append(
                {
                    "p_key": point["p_key"],
                    "p_value": point["p_value"],
                    "driver_class": driver,
                    "runs": len(seeds),
                    "seed_min": min(seeds),
                    "seed_max": max(seeds),
                    "seed_set_sha256": hashlib.sha256(encoded).hexdigest(),
                    "seed_source": source,
                    "command_p_fully_checked": p_checked,
                }
            )
    reference_key = (points[0]["p_key"], DRIVERS[0])
    reference = seed_sets[reference_key]
    mismatches = [key for key, seeds in seed_sets.items() if seeds != reference]
    if mismatches:
        raise ValueError(
            f"Sweep arms do not use one matched seed set; reference={reference_key}, "
            f"mismatches={mismatches}"
        )
    return seed_sets, pd.DataFrame(summary_rows)


def weak_connectivity(edges: Iterable[tuple[str, str]]) -> tuple[int, int, int]:
    parent: dict[str, str] = {}
    size: dict[str, int] = {}

    def find(node: str) -> str:
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != node:
            following = parent[node]
            parent[node] = root
            node = following
        return root

    def add(node: str) -> None:
        if node not in parent:
            parent[node] = node
            size[node] = 1

    for source, target in edges:
        add(source)
        add(target)
        left = find(source)
        right = find(target)
        if left == right:
            continue
        if size[left] < size[right]:
            left, right = right, left
        parent[right] = left
        size[left] += size[right]
    if not parent:
        return 0, 0, 0
    roots = [node for node in parent if parent[node] == node]
    return len(parent), len(roots), max(size[root] for root in roots)


def distribution_metrics(values: pd.Series, prefix: str) -> dict[str, float]:
    if values.empty:
        return {
            f"{prefix}_mean": math.nan,
            f"{prefix}_median": math.nan,
            f"{prefix}_q25": math.nan,
            f"{prefix}_q75": math.nan,
            f"{prefix}_min": math.nan,
            f"{prefix}_max": math.nan,
        }
    return {
        f"{prefix}_mean": float(values.mean()),
        f"{prefix}_median": float(values.median()),
        f"{prefix}_q25": float(values.quantile(0.25)),
        f"{prefix}_q75": float(values.quantile(0.75)),
        f"{prefix}_min": float(values.min()),
        f"{prefix}_max": float(values.max()),
    }


def summarize_network(
    frame: pd.DataFrame,
    candidates: list[str],
    expression_ids: set[str],
    *,
    prefix: dict[str, Any],
    seed_runs: int,
) -> tuple[dict[str, Any], pd.Series]:
    source_set = set(frame["source"])
    unexpected = source_set - set(candidates)
    if unexpected:
        preview = sorted(unexpected)[:5]
        raise ValueError(f"Non-candidate sources for {prefix}: {preview}")
    incident_set = source_set | set(frame["target"])
    unexpected_incident = incident_set - expression_ids
    if unexpected_incident:
        preview = sorted(unexpected_incident)[:5]
        raise ValueError(
            f"Edges contain IDs absent from the expression matrix for {prefix}: "
            f"{preview}"
        )
    if not set(candidates).issubset(expression_ids):
        raise ValueError(
            f"Candidate list contains IDs absent from the expression matrix for {prefix}"
        )
    if seed_runs > 0 and len(frame):
        if (frame["support_count"] > seed_runs).any() or not np.allclose(
            frame["support_fraction"],
            frame["support_count"] / seed_runs,
            rtol=0.0,
            atol=1e-12,
        ):
            raise ValueError(
                f"Support counts/fractions disagree with {seed_runs} matched seeds for {prefix}"
            )
    target_sizes = (
        frame.groupby(level="source", sort=False)
        .size()
        .reindex(candidates, fill_value=0)
    )
    incident_nodes, components, largest = weak_connectivity(frame.index)
    row: dict[str, Any] = {
        **prefix,
        "seed_runs": seed_runs,
        "candidate_drivers": len(candidates),
        "active_drivers": int((target_sizes > 0).sum()),
        "active_driver_fraction": float((target_sizes > 0).mean()),
        "consensus_edges": len(frame),
        "incident_nodes": incident_nodes,
        "expression_nodes": len(expression_ids),
        "incident_node_fraction_expression": incident_nodes / len(expression_ids),
        "weak_components": components,
        "largest_weak_component_nodes": largest,
        "largest_weak_component_fraction_incident": (
            largest / incident_nodes if incident_nodes else 0.0
        ),
        "target_size_mean_zero_filled": float(target_sizes.mean()),
        "target_size_median_zero_filled": float(target_sizes.median()),
        "target_size_q25_zero_filled": float(target_sizes.quantile(0.25)),
        "target_size_q75_zero_filled": float(target_sizes.quantile(0.75)),
        "target_size_max": int(target_sizes.max()),
        **distribution_metrics(frame["support_fraction"], "support_fraction"),
        **distribution_metrics(frame["support_count"], "support_count"),
        **distribution_metrics(frame["MI"], "consensus_MI"),
    }
    return row, target_sizes


def safe_ratio(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else math.nan


def safe_correlation(left: pd.Series, right: pd.Series, method: str) -> float:
    if len(left) < 3 or left.nunique() < 2 or right.nunique() < 2:
        return math.nan
    return float(left.corr(right, method=method))


def overlap_row(
    reference: pd.DataFrame,
    comparison: pd.DataFrame,
    reference_sizes: pd.Series,
    comparison_sizes: pd.Series,
    *,
    prefix: dict[str, Any],
) -> dict[str, Any]:
    reference_edges = set(reference.index)
    comparison_edges = set(comparison.index)
    common = reference_edges & comparison_edges
    union = reference_edges | comparison_edges
    reference_active = set(reference["source"])
    comparison_active = set(comparison["source"])
    active_union = reference_active | comparison_active
    common_index = pd.MultiIndex.from_tuples(sorted(common), names=["source", "target"])
    if common:
        reference_mi = reference.loc[common_index, "MI"].reset_index(drop=True)
        comparison_mi = comparison.loc[common_index, "MI"].reset_index(drop=True)
        reference_support = reference.loc[common_index, "support_fraction"].reset_index(drop=True)
        comparison_support = comparison.loc[common_index, "support_fraction"].reset_index(drop=True)
    else:
        reference_mi = pd.Series(dtype=float)
        comparison_mi = pd.Series(dtype=float)
        reference_support = pd.Series(dtype=float)
        comparison_support = pd.Series(dtype=float)
    return {
        **prefix,
        "reference_edges": len(reference_edges),
        "comparison_edges": len(comparison_edges),
        "intersection_edges": len(common),
        "union_edges": len(union),
        "edge_jaccard": safe_ratio(len(common), len(union)),
        "reference_edge_retention": safe_ratio(len(common), len(reference_edges)),
        "comparison_edge_shared_fraction": safe_ratio(len(common), len(comparison_edges)),
        "lost_from_reference": len(reference_edges - comparison_edges),
        "gained_in_comparison": len(comparison_edges - reference_edges),
        "reference_active_drivers": len(reference_active),
        "comparison_active_drivers": len(comparison_active),
        "active_driver_intersection": len(reference_active & comparison_active),
        "active_driver_union": len(active_union),
        "active_driver_jaccard": safe_ratio(
            len(reference_active & comparison_active), len(active_union)
        ),
        "target_size_pearson_all": safe_correlation(reference_sizes, comparison_sizes, "pearson"),
        "target_size_spearman_all": safe_correlation(reference_sizes, comparison_sizes, "spearman"),
        "common_edge_mi_n": len(common),
        "common_edge_mi_pearson": safe_correlation(reference_mi, comparison_mi, "pearson"),
        "common_edge_mi_spearman": safe_correlation(reference_mi, comparison_mi, "spearman"),
        "common_edge_support_pearson": safe_correlation(
            reference_support, comparison_support, "pearson"
        ),
        "common_edge_support_spearman": safe_correlation(
            reference_support, comparison_support, "spearman"
        ),
    }


def locate_anchor(points: list[dict[str, Any]], requested: float) -> dict[str, Any]:
    matches = [
        point
        for point in points
        if math.isclose(point["p_value"], requested, rel_tol=1e-12, abs_tol=0.0)
    ]
    if len(matches) != 1:
        available = ", ".join(f"{point['p_value']:.12g}" for point in points)
        raise ValueError(f"Anchor p={requested:g} is absent or ambiguous; available: {available}")
    return matches[0]


def locate_pr66_stage(work_root: Path) -> Path:
    results_root = work_root / "results"
    candidates: list[Path] = []
    if all((results_root / driver).is_dir() for driver in DRIVERS):
        candidates.append(results_root)
    if results_root.is_dir():
        candidates.extend(
            path
            for path in sorted(results_root.iterdir())
            if path.is_dir()
            and path.name.lower().startswith("pr66")
            and all((path / driver).is_dir() for driver in DRIVERS)
        )
    unique = list(dict.fromkeys(candidates))
    if len(unique) != 1:
        raise ValueError(
            f"Could not identify exactly one PR66 results stage under {work_root}; "
            f"candidates={unique}"
        )
    return unique[0]


def save_plot(figure: plt.Figure, plots_root: Path, stem: str) -> None:
    figure.savefig(plots_root / f"{stem}.png", dpi=180)
    figure.savefig(plots_root / f"{stem}.svg", metadata=SVG_METADATA)


def plot_core_metrics(
    summary: pd.DataFrame,
    points: list[dict[str, Any]],
    anchor: dict[str, Any],
    plots_root: Path,
    context_summary: pd.DataFrame | None,
) -> None:
    metrics = (
        ("consensus_edges", "Consensus edges", "symlog"),
        ("active_driver_fraction", "Active candidate drivers", "fraction"),
        (
            "largest_weak_component_fraction_incident",
            "Largest weak component / incident nodes",
            "fraction",
        ),
        ("support_fraction_median", "Median consensus support", "fraction"),
        ("consensus_MI_median", "Median consensus MI", "linear"),
    )
    x_values = np.asarray([point["log10_p"] for point in points])
    x_labels = [point["p_label"] for point in points]
    figure, axes = plt.subplots(
        len(DRIVERS), len(metrics), figsize=(19.0, 7.2), constrained_layout=True, squeeze=False
    )
    colors = {"tf": "#0072B2", "sig": "#D55E00"}
    for row_index, driver in enumerate(DRIVERS):
        subset = summary[summary["driver_class"] == driver].sort_values("p_value")
        if len(subset) != len(points):
            raise ValueError(f"Incomplete summary rows for {driver}")
        for column_index, (metric, label, scale) in enumerate(metrics):
            axis = axes[row_index, column_index]
            y_values = subset[metric].to_numpy(dtype=float)
            axis.plot(
                x_values,
                y_values,
                marker="o",
                markersize=4.5,
                linewidth=1.8,
                color=colors[driver],
                label="PR67 sweep",
            )
            if metric in ("support_fraction_median", "consensus_MI_median"):
                prefix = "support_fraction" if metric.startswith("support") else "consensus_MI"
                axis.fill_between(
                    x_values,
                    subset[f"{prefix}_q25"].to_numpy(dtype=float),
                    subset[f"{prefix}_q75"].to_numpy(dtype=float),
                    color=colors[driver],
                    alpha=0.14,
                    linewidth=0.0,
                    label="IQR",
                )
            if context_summary is not None:
                context_value = float(
                    context_summary.loc[context_summary["driver_class"] == driver, metric].iloc[0]
                )
                axis.axhline(
                    context_value,
                    color="#555555",
                    linestyle="--",
                    linewidth=1.2,
                    label="PR66 context",
                )
            axis.axvline(anchor["log10_p"], color="#777777", linestyle=":", linewidth=1.0)
            if scale == "symlog":
                axis.set_yscale("symlog", linthresh=1)
            elif scale == "fraction":
                axis.set_ylim(-0.025, 1.025)
            axis.set_title(label)
            axis.grid(alpha=0.22)
            if row_index == len(DRIVERS) - 1:
                axis.set_xticks(x_values, x_labels, rotation=35, ha="right")
                axis.set_xlabel("Per-subsample p (x position = log10(p))")
            else:
                axis.set_xticks(x_values, [])
            if column_index == 0:
                axis.set_ylabel(driver.upper())
            if row_index == 0 and column_index == 0:
                axis.legend(frameon=False, fontsize=8)
    figure.suptitle(
        "BRCA100 PR67 per-subsample threshold sweep; dotted line is the "
        f"p={anchor['p_value']:.3g} anchor"
    )
    save_plot(figure, plots_root, "core_metrics_vs_log10_p")
    plt.close(figure)


def plot_overlaps(
    adjacent: pd.DataFrame,
    anchor_overlap: pd.DataFrame,
    plots_root: Path,
    anchor_p: float,
) -> None:
    figure, axes = plt.subplots(1, len(DRIVERS), figsize=(12.5, 4.8), constrained_layout=True)
    for axis, driver in zip(axes, DRIVERS):
        anchored = anchor_overlap[anchor_overlap["driver_class"] == driver].sort_values(
            "comparison_p_value"
        )
        neighboring = adjacent[adjacent["driver_class"] == driver].sort_values(
            "comparison_p_value"
        )
        axis.plot(
            anchored["comparison_log10_p"],
            anchored["edge_jaccard"],
            marker="o",
            color="#0072B2",
            label=f"vs p={anchor_p:.3g} anchor",
        )
        axis.plot(
            neighboring["comparison_log10_p"],
            neighboring["edge_jaccard"],
            marker="s",
            color="#D55E00",
            label="vs adjacent stricter p",
        )
        axis.set_ylim(-0.025, 1.025)
        axis.set_xlabel("log10(per-subsample p)")
        axis.set_ylabel("Directed-edge Jaccard")
        axis.set_title(driver.upper())
        axis.grid(alpha=0.22)
    axes[0].legend(frameon=False)
    figure.suptitle("Consensus-network edge stability across the PR67 threshold sweep")
    save_plot(figure, plots_root, "edge_overlap_vs_log10_p")
    plt.close(figure)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, float_format="%.12g")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-root",
        type=Path,
        required=True,
        help="Sweep root containing inputs/ and results/<p_key>/point_manifest.json",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Output directory (default: <work-root>/results/analysis)",
    )
    parser.add_argument(
        "--anchor-p",
        type=float,
        default=1e-7,
        help="Sweep p used for anchor-overlap comparisons (default: 1e-7)",
    )
    parser.add_argument(
        "--pr66-work-root",
        type=Path,
        help="Optional prior matched BRCA100 work root used only as plot/table context",
    )
    args = parser.parse_args()

    work_root = args.work_root.resolve()
    results_root = work_root / "results"
    if not results_root.is_dir():
        raise ValueError(f"Missing results directory: {results_root}")
    output_root = (
        args.output_root.resolve()
        if args.output_root is not None
        else results_root / "analysis"
    )
    plots_root = output_root / "plots"
    output_root.mkdir(parents=True, exist_ok=True)
    plots_root.mkdir(parents=True, exist_ok=True)

    points = point_metadata(results_root)
    anchor = locate_anchor(points, args.anchor_p)
    candidates = candidate_lists(work_root / "inputs")
    expression_path = work_root / "inputs" / "BRCA100.exp"
    expression_ids = expression_identifiers(expression_path)
    seed_sets, seed_summary = collect_seed_sets(work_root, points)
    write_tsv(seed_summary, output_root / "seed_manifest_summary.tsv")

    point_table = pd.DataFrame(
        [
            {
                key: point[key]
                for key in (
                    "p_key",
                    "p_value",
                    "log10_p",
                    "p_token",
                    "p_label",
                    "role",
                    "mi_cutoff",
                    "validation_class",
                    "tail_extrapolated",
                    "commit",
                    "model_sha256",
                    "manifest_sha256",
                )
            }
            for point in points
        ]
    )
    write_tsv(point_table, output_root / "point_manifest_summary.tsv")

    networks: dict[tuple[str, str], pd.DataFrame] = {}
    sizes: dict[tuple[str, str], pd.Series] = {}
    summary_rows: list[dict[str, Any]] = []
    for point in points:
        for driver in DRIVERS:
            key = (point["p_key"], driver)
            frame = merge_network_support(point["root"] / driver)
            networks[key] = frame
            row, target_sizes = summarize_network(
                frame,
                candidates[driver],
                expression_ids,
                prefix={
                    "p_key": point["p_key"],
                    "p_value": point["p_value"],
                    "log10_p": point["log10_p"],
                    "p_token": point["p_token"],
                    "p_label": point["p_label"],
                    "role": point["role"],
                    "mi_cutoff": point["mi_cutoff"],
                    "validation_class": point["validation_class"],
                    "tail_extrapolated": point["tail_extrapolated"],
                    "driver_class": driver,
                },
                seed_runs=len(seed_sets[key]),
            )
            summary_rows.append(row)
            sizes[key] = target_sizes
    summary = pd.DataFrame(summary_rows).sort_values(["driver_class", "p_value"])
    write_tsv(summary, output_root / "network_summary.tsv")

    adjacent_rows: list[dict[str, Any]] = []
    for reference_point, comparison_point in zip(points[:-1], points[1:]):
        for driver in DRIVERS:
            reference_key = (reference_point["p_key"], driver)
            comparison_key = (comparison_point["p_key"], driver)
            adjacent_rows.append(
                overlap_row(
                    networks[reference_key],
                    networks[comparison_key],
                    sizes[reference_key],
                    sizes[comparison_key],
                    prefix={
                        "driver_class": driver,
                        "reference_p_key": reference_point["p_key"],
                        "reference_p_value": reference_point["p_value"],
                        "reference_log10_p": reference_point["log10_p"],
                        "comparison_p_key": comparison_point["p_key"],
                        "comparison_p_value": comparison_point["p_value"],
                        "comparison_log10_p": comparison_point["log10_p"],
                    },
                )
            )
    adjacent = pd.DataFrame(adjacent_rows)
    write_tsv(adjacent, output_root / "adjacent_overlap.tsv")

    anchor_rows: list[dict[str, Any]] = []
    for comparison_point in points:
        for driver in DRIVERS:
            reference_key = (anchor["p_key"], driver)
            comparison_key = (comparison_point["p_key"], driver)
            anchor_rows.append(
                overlap_row(
                    networks[reference_key],
                    networks[comparison_key],
                    sizes[reference_key],
                    sizes[comparison_key],
                    prefix={
                        "driver_class": driver,
                        "reference_p_key": anchor["p_key"],
                        "reference_p_value": anchor["p_value"],
                        "reference_log10_p": anchor["log10_p"],
                        "comparison_p_key": comparison_point["p_key"],
                        "comparison_p_value": comparison_point["p_value"],
                        "comparison_log10_p": comparison_point["log10_p"],
                    },
                )
            )
    anchor_overlap = pd.DataFrame(anchor_rows)
    write_tsv(anchor_overlap, output_root / "anchor_overlap.tsv")

    context_summary: pd.DataFrame | None = None
    context_provenance: dict[str, Any] | None = None
    if args.pr66_work_root is not None:
        pr66_work_root = args.pr66_work_root.resolve()
        pr66_stage = locate_pr66_stage(pr66_work_root)
        prior_expression_path = pr66_work_root / "inputs" / "BRCA100.exp"
        if prior_expression_path.is_file() and sha256_file(prior_expression_path) != sha256_file(
            expression_path
        ):
            raise ValueError("PR66 and sweep expression matrices differ")
        context_networks: dict[str, pd.DataFrame] = {}
        context_sizes: dict[str, pd.Series] = {}
        context_rows: list[dict[str, Any]] = []
        for driver in DRIVERS:
            prior_candidates_path = pr66_work_root / "inputs" / CANDIDATE_FILES[driver]
            current_candidates_path = work_root / "inputs" / CANDIDATE_FILES[driver]
            if (
                prior_candidates_path.is_file()
                and sha256_file(prior_candidates_path)
                != sha256_file(current_candidates_path)
            ):
                raise ValueError(f"PR66 and sweep candidate lists differ for {driver}")
            frame = merge_network_support(pr66_stage / driver)
            context_networks[driver] = frame
            row, target_sizes = summarize_network(
                frame,
                candidates[driver],
                expression_ids,
                prefix={"context": "pr66", "driver_class": driver},
                seed_runs=0,
            )
            context_rows.append(row)
            context_sizes[driver] = target_sizes
        context_summary = pd.DataFrame(context_rows)
        write_tsv(context_summary, output_root / "pr66_context_summary.tsv")
        context_overlap_rows: list[dict[str, Any]] = []
        for comparison_point in points:
            for driver in DRIVERS:
                comparison_key = (comparison_point["p_key"], driver)
                context_overlap_rows.append(
                    overlap_row(
                        context_networks[driver],
                        networks[comparison_key],
                        context_sizes[driver],
                        sizes[comparison_key],
                        prefix={
                            "driver_class": driver,
                            "reference_context": "pr66",
                            "comparison_p_key": comparison_point["p_key"],
                            "comparison_p_value": comparison_point["p_value"],
                            "comparison_log10_p": comparison_point["log10_p"],
                        },
                    )
                )
        write_tsv(pd.DataFrame(context_overlap_rows), output_root / "pr66_context_overlap.tsv")
        context_provenance = {
            "work_root": str(pr66_work_root),
            "stage_root": str(pr66_stage),
        }

    plot_core_metrics(summary, points, anchor, plots_root, context_summary)
    plot_overlaps(adjacent, anchor_overlap, plots_root, anchor["p_value"])

    analysis_manifest = {
        "schema": "sjaracne-brca100-pr67-threshold-sweep-analysis-v1",
        "work_root": str(work_root),
        "output_root": str(output_root),
        "anchor_p": anchor["p_value"],
        "anchor_p_key": anchor["p_key"],
        "p_keys_in_increasing_p_order": [point["p_key"] for point in points],
        "candidate_files": {
            driver: {
                "path": str(work_root / "inputs" / filename),
                "sha256": sha256_file(work_root / "inputs" / filename),
                "count": len(candidates[driver]),
            }
            for driver, filename in CANDIDATE_FILES.items()
        },
        "expression_file": {
            "path": str(expression_path),
            "sha256": sha256_file(expression_path),
            "row_count": len(expression_ids),
        },
        "matched_seed_count": len(next(iter(seed_sets.values()))),
        "pr66_context": context_provenance,
        "interpretation_scope": (
            "Network/statistical QC only; no biological-reference or downstream "
            "NetBID reproducibility criterion was applied."
        ),
    }
    (output_root / "analysis_manifest.json").write_text(
        json.dumps(analysis_manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"Wrote threshold-sweep analysis to {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
