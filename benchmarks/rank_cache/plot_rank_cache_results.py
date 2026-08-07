#!/usr/bin/env python3
"""Create publication-ready plots for the AP-MI rank-cache benchmark."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


HERE = Path(__file__).resolve().parent
DEFAULT_SUMMARY = HERE / "results_2026-08-07.csv"
DEFAULT_RAW = HERE / "raw_timings_2026-08-07.csv"
DEFAULT_OUTPUT_DIR = HERE / "figures"

LEGACY_COLOR = "#D55E00"
CACHE_COLOR = "#0072B2"
MEMORY_COLOR = "#7A5195"
TEXT_COLOR = "#24292F"
GRID_COLOR = "#D8DEE4"

SWEEPS = (
    {
        "title": "Hub-count sweep",
        "fixed": "G = 5,000; N = 1,000",
        "x_label": "Hub genes (H)",
        "x_values": (1, 10, 50, 100),
        "cases": (
            "h1_g5000_n1000",
            "h10_g5000_n1000",
            "h50_g5000_n1000",
            "h100_g5000_n1000",
        ),
        "labels": ("1", "10", "50", "100"),
    },
    {
        "title": "Observation-count sweep",
        "fixed": "G = 5,000; H = 100",
        "x_label": "Observations (N)",
        "x_values": (100, 250, 500, 1000),
        "cases": (
            "h100_g5000_n0100",
            "h100_g5000_n0250",
            "h100_g5000_n0500",
            "h100_g5000_n1000",
        ),
        "labels": ("100", "250", "500", "1,000"),
    },
    {
        "title": "Gene-count sweep",
        "fixed": "N = 500; H = 100",
        "x_label": "Genes (G)",
        "x_values": (1000, 5000, 10000, 19936),
        "cases": (
            "h100_g1000_n500",
            "h100_g5000_n0500",
            "h100_g10000_n500",
            "h100_g19936_n500",
        ),
        "labels": ("1,000", "5,000", "10,000", "19,936†"),
    },
)


@dataclass
class RawMeasurements:
    repeat_ids: list[int] = field(default_factory=list)
    wall_s: list[float] = field(default_factory=list)
    max_rss_kb: list[float] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY,
        help="Summary benchmark CSV (default: %(default)s)",
    )
    parser.add_argument(
        "--raw",
        type=Path,
        default=DEFAULT_RAW,
        help="Raw timing CSV (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Figure output directory (default: %(default)s)",
    )
    return parser.parse_args()


def read_summary(path: Path) -> dict[str, dict[str, float]]:
    required = {
        "case",
        "repetitions",
        "baseline_wall_median_s",
        "optimized_wall_median_s",
        "wall_speedup_x",
        "baseline_rss_median_kb",
        "optimized_rss_median_kb",
    }
    rows: dict[str, dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        for raw_row in reader:
            case = raw_row["case"]
            if case in rows:
                raise ValueError(f"Duplicate summary case: {case}")
            rows[case] = {
                key: float(raw_row[key]) for key in required if key != "case"
            }
            if not all(math.isfinite(value) for value in rows[case].values()):
                raise ValueError(f"Non-finite summary value for case: {case}")
    return rows


def read_raw(path: Path) -> dict[tuple[str, str], RawMeasurements]:
    required = {
        "implementation",
        "case",
        "repeat",
        "wall_s",
        "max_rss_kb",
        "exit_status",
    }
    rows: dict[tuple[str, str], RawMeasurements] = defaultdict(RawMeasurements)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        for row in reader:
            implementation = row["implementation"]
            if implementation not in {"baseline", "optimized"}:
                raise ValueError(f"Unexpected implementation: {implementation}")
            if int(row["exit_status"]) != 0:
                raise ValueError(
                    f"Nonzero benchmark exit for {implementation}/{row['case']}"
                )
            repeat_id = int(row["repeat"])
            wall_s = float(row["wall_s"])
            max_rss_kb = float(row["max_rss_kb"])
            if not math.isfinite(wall_s) or wall_s <= 0:
                raise ValueError(
                    f"Invalid wall time for {implementation}/{row['case']}: {wall_s}"
                )
            if not math.isfinite(max_rss_kb) or max_rss_kb <= 0:
                raise ValueError(
                    f"Invalid peak RSS for {implementation}/{row['case']}: "
                    f"{max_rss_kb}"
                )
            measurements = rows[(implementation, row["case"])]
            if repeat_id in measurements.repeat_ids:
                raise ValueError(
                    f"Duplicate repeat {repeat_id} for {implementation}/{row['case']}"
                )
            measurements.repeat_ids.append(repeat_id)
            measurements.wall_s.append(wall_s)
            measurements.max_rss_kb.append(max_rss_kb)
    return dict(rows)


def validate_inputs(
    summary: dict[str, dict[str, float]],
    raw: dict[tuple[str, str], RawMeasurements],
) -> None:
    cases = {case for sweep in SWEEPS for case in sweep["cases"]}
    missing_summary = cases.difference(summary)
    if missing_summary:
        raise ValueError(f"Summary is missing cases: {sorted(missing_summary)}")
    unexpected_summary = set(summary).difference(cases)
    if unexpected_summary:
        raise ValueError(f"Summary has unexpected cases: {sorted(unexpected_summary)}")
    unexpected_raw = {case for _, case in raw}.difference(cases)
    if unexpected_raw:
        raise ValueError(f"Raw timings have unexpected cases: {sorted(unexpected_raw)}")

    for case in sorted(cases):
        expected_repetitions = int(summary[case]["repetitions"])
        for implementation, summary_key in (
            ("baseline", "baseline_wall_median_s"),
            ("optimized", "optimized_wall_median_s"),
        ):
            measurements = raw.get((implementation, case))
            if not measurements:
                raise ValueError(f"Raw timings missing for {implementation}/{case}")
            if sorted(measurements.repeat_ids) != list(
                range(1, expected_repetitions + 1)
            ):
                raise ValueError(
                    f"Unexpected repeats for {implementation}/{case}: "
                    f"{sorted(measurements.repeat_ids)}"
                )
            observed = median(measurements.wall_s)
            expected = summary[case][summary_key]
            if not math.isclose(observed, expected, abs_tol=0.005):
                raise ValueError(
                    f"Median mismatch for {implementation}/{case}: "
                    f"raw={observed}, summary={expected}"
                )

        baseline = raw[("baseline", case)]
        optimized = raw[("optimized", case)]
        observed_speedup = median(baseline.wall_s) / median(optimized.wall_s)
        if not math.isclose(
            observed_speedup, summary[case]["wall_speedup_x"], abs_tol=0.01
        ):
            raise ValueError(f"Speedup mismatch for case: {case}")
        for measurements, summary_key in (
            (baseline, "baseline_rss_median_kb"),
            (optimized, "optimized_rss_median_kb"),
        ):
            if not math.isclose(
                median(measurements.max_rss_kb),
                summary[case][summary_key],
                abs_tol=0.5,
            ):
                raise ValueError(f"Peak-RSS median mismatch for case: {case}")


def set_numeric_x_axis(ax: plt.Axes, sweep: dict[str, object]) -> None:
    x_values = sweep["x_values"]
    ax.set_xscale("log")
    ax.set_xlim(x_values[0] / 1.18, x_values[-1] * 1.18)
    ax.set_xticks(x_values, sweep["labels"])
    ax.minorticks_off()
    ax.set_xlabel(sweep["x_label"], labelpad=8)


def configure_style() -> None:
    matplotlib.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#8C959F",
            "axes.labelcolor": TEXT_COLOR,
            "axes.titlecolor": TEXT_COLOR,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "text.color": TEXT_COLOR,
            "xtick.color": TEXT_COLOR,
            "ytick.color": TEXT_COLOR,
            "legend.fontsize": 9.5,
            "svg.fonttype": "none",
            "svg.hashsalt": "sjaracne-rank-cache-v1",
        }
    )


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.8, alpha=0.75)
    ax.set_axisbelow(True)


def add_panel_header(ax: plt.Axes, title: str, fixed: str) -> None:
    ax.set_title(title, pad=20)
    ax.text(
        0.5,
        1.01,
        fixed,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#57606A",
    )


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {"Creator": "SJARACNe rank-cache benchmark plotting script"}
    svg_path = output_dir / f"{stem}.svg"
    fig.savefig(
        output_dir / f"{stem}.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": metadata["Creator"]},
    )
    fig.savefig(
        svg_path,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Creator": metadata["Creator"], "Date": None},
    )
    # Matplotlib emits trailing spaces in multiline SVG path data. Normalize the
    # generated text so repository whitespace checks remain useful.
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    plt.close(fig)


def plot_runtime(
    raw: dict[tuple[str, str], RawMeasurements],
    output_dir: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 5.05))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.68, bottom=0.22, wspace=0.27)
    fig.suptitle(
        "Bootstrap rank caching reduces AP–MI runtime",
        x=0.5,
        y=0.965,
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.905,
        "Lines show median wall time; small dots show individual runs",
        ha="center",
        fontsize=9.5,
        color="#57606A",
    )

    for ax, sweep in zip(axes, SWEEPS):
        cases = sweep["cases"]
        x = sweep["x_values"]
        legacy = [median(raw[("baseline", case)].wall_s) for case in cases]
        cached = [median(raw[("optimized", case)].wall_s) for case in cases]

        ax.plot(
            x,
            legacy,
            color=LEGACY_COLOR,
            linewidth=2.3,
            marker="o",
            markersize=6.5,
            label="Legacy",
            zorder=3,
        )
        ax.plot(
            x,
            cached,
            color=CACHE_COLOR,
            linewidth=2.3,
            marker="o",
            markersize=6.5,
            label="Rank cache",
            zorder=3,
        )

        for x_value, case in zip(x, cases):
            for implementation, color, center in (
                ("baseline", LEGACY_COLOR, x_value * 0.965),
                ("optimized", CACHE_COLOR, x_value * 1.035),
            ):
                values = raw[(implementation, case)].wall_s
                if len(values) == 1:
                    offsets = [0.0]
                else:
                    offsets = [
                        -0.008 + 0.016 * index / (len(values) - 1)
                        for index in range(len(values))
                    ]
                ax.scatter(
                    [center * (1 + offset) for offset in offsets],
                    values,
                    s=18,
                    color=color,
                    alpha=0.42,
                    edgecolors="none",
                    zorder=2,
                )

        speedup = legacy[-1] / cached[-1]
        ax.annotate(
            f"{speedup:.2f}× faster",
            (x[-1], cached[-1]),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            color=CACHE_COLOR,
            fontsize=9,
            fontweight="bold",
        )
        set_numeric_x_axis(ax, sweep)
        ax.set_ylim(0, max(legacy) * 1.16)
        add_panel_header(ax, sweep["title"], sweep["fixed"])
        style_axis(ax)

    axes[0].set_ylabel("Wall time (seconds)")
    legend_handles = [
        Line2D([0], [0], color=LEGACY_COLOR, marker="o", lw=2.3, label="Legacy"),
        Line2D([0], [0], color=CACHE_COLOR, marker="o", lw=2.3, label="Rank cache"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.835),
        ncol=2,
        frameon=False,
    )
    fig.text(
        0.5,
        0.055,
        "Medians of 3 sequential runs; †G = 19,936 is one stress run. "
        "Log-scaled x axes. AP–MI kernel benchmark, not the complete workflow.",
        ha="center",
        fontsize=8.3,
        color="#57606A",
    )
    save_figure(fig, output_dir, "rank_cache_runtime_sweeps")


def plot_speedup(
    raw: dict[tuple[str, str], RawMeasurements], output_dir: Path
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.6), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.76, bottom=0.23, wspace=0.17)
    fig.suptitle(
        "Measured AP–MI speedup from bootstrap rank caching",
        x=0.5,
        y=0.965,
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.905,
        "Baseline median ÷ rank-cache median; higher is faster",
        ha="center",
        fontsize=9.5,
        color="#57606A",
    )

    for ax, sweep in zip(axes, SWEEPS):
        cases = sweep["cases"]
        x = sweep["x_values"]
        values = [
            median(raw[("baseline", case)].wall_s)
            / median(raw[("optimized", case)].wall_s)
            for case in cases
        ]
        ax.axhline(1.0, color="#8C959F", linestyle="--", linewidth=1.2, zorder=1)
        ax.plot(
            x,
            values,
            color=CACHE_COLOR,
            linewidth=2.6,
            marker="o",
            markersize=7.2,
            zorder=3,
        )
        for position, value in zip(x, values):
            ax.annotate(
                f"{value:.2f}×",
                (position, value),
                xytext=(0, 9),
                textcoords="offset points",
                ha="center",
                color=CACHE_COLOR,
                fontsize=9,
                fontweight="bold",
            )
        set_numeric_x_axis(ax, sweep)
        ax.set_ylim(0.75, 5.75)
        add_panel_header(ax, sweep["title"], sweep["fixed"])
        style_axis(ax)

    axes[0].set_ylabel("Speedup (fold)")
    fig.text(
        0.5,
        0.055,
        "Medians of 3 sequential runs; †G = 19,936 is one stress run. "
        "Log-scaled x axes; dashed line marks no improvement (1×).",
        ha="center",
        fontsize=8.3,
        color="#57606A",
    )
    save_figure(fig, output_dir, "rank_cache_speedup_sweeps")


def plot_memory(
    raw: dict[tuple[str, str], RawMeasurements], output_dir: Path
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.5), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.76, bottom=0.23, wspace=0.17)
    fig.suptitle(
        "Memory cost of bootstrap rank caching",
        x=0.5,
        y=0.965,
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.905,
        "Increase in median peak resident memory relative to the legacy implementation",
        ha="center",
        fontsize=9.5,
        color="#57606A",
    )

    for ax, sweep in zip(axes, SWEEPS):
        cases = sweep["cases"]
        x = sweep["x_values"]
        values = [
            100
            * (
                median(raw[("optimized", case)].max_rss_kb)
                / median(raw[("baseline", case)].max_rss_kb)
                - 1
            )
            for case in cases
        ]
        ax.vlines(x, 0, values, color=MEMORY_COLOR, linewidth=6, alpha=0.38)
        ax.plot(
            x,
            values,
            color=MEMORY_COLOR,
            linewidth=2.1,
            marker="o",
            markersize=7,
            zorder=3,
        )
        for index, (x_value, value) in enumerate(zip(x, values)):
            ax.annotate(
                f"+{value:.1f}%",
                (x_value, value),
                xytext=(0, 5 + 10 * (index % 2)),
                textcoords="offset points",
                ha="center",
                color=MEMORY_COLOR,
                fontsize=8.7,
                fontweight="bold",
            )
        set_numeric_x_axis(ax, sweep)
        ax.set_ylim(0, 21.5)
        add_panel_header(ax, sweep["title"], sweep["fixed"])
        style_axis(ax)

    axes[0].set_ylabel("Peak RSS increase (%)")
    fig.text(
        0.5,
        0.055,
        "Medians of 3 sequential runs; †G = 19,936 is one stress run. "
        "Log-scaled x axes; cache stores one 32-bit rank per cached gene-observation pair.",
        ha="center",
        fontsize=8.3,
        color="#57606A",
    )
    save_figure(fig, output_dir, "rank_cache_memory_overhead")


def main() -> None:
    args = parse_args()
    summary = read_summary(args.summary)
    raw = read_raw(args.raw)
    validate_inputs(summary, raw)
    configure_style()
    plot_runtime(raw, args.output_dir)
    plot_speedup(raw, args.output_dir)
    plot_memory(raw, args.output_dir)
    print(f"Wrote rank-cache benchmark figures to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
