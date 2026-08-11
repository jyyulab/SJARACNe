#!/usr/bin/env python3

import argparse
import csv
import hashlib
import statistics
import subprocess
import time
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark removal of SJARACNe's unused bandwidth pass."
    )
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--rank-data", type=Path, required=True)
    parser.add_argument("--brca-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    return parser.parse_args()


def network_data_hash(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for line in handle:
            if line.strip() and not line.startswith(b">"):
                digest.update(line)
    return digest.hexdigest()


def cases(rank_data, brca_data, empty_adjacency):
    common = ["-p", "1", "-e", "0"]
    return {
        "preprocess_g01000_n0500": [
            "-i", str(rank_data / "expression_g01000_n0500.exp"),
            "-j", str(empty_adjacency), *common,
        ],
        "preprocess_g05000_n1000": [
            "-i", str(rank_data / "expression_g05000_n1000.exp"),
            "-j", str(empty_adjacency), *common,
        ],
        "preprocess_g19936_n0500": [
            "-i", str(rank_data / "expression_g19936_n0500.exp"),
            "-j", str(empty_adjacency), *common,
        ],
        "ap_h0001_g05000_n1000": [
            "-i", str(rank_data / "expression_g05000_n1000.exp"),
            "-s", str(rank_data / "hubs_h0001.txt"),
            "-S", "17", "-r", "1", "-t", "100", "-e", "1",
        ],
        "ap_h0100_g05000_n1000": [
            "-i", str(rank_data / "expression_g05000_n1000.exp"),
            "-s", str(rank_data / "hubs_h0100.txt"),
            "-S", "17", "-r", "1", "-t", "100", "-e", "1",
        ],
        "brca100_adjacency_replay": [
            "-i", str(brca_data / "BRCA100.exp"),
            "-j", str(brca_data / "TF_run.adj"),
            "-s", str(brca_data / "BRCA100_TF.txt"),
            "-l", str(brca_data / "BRCA100_TF.txt"),
            "-p", "1", "-e", "0",
        ],
    }


def run_once(binary, arguments, output, metrics, stdout, stderr):
    command = [
        "/usr/bin/time", "-f", "%e\t%U\t%S\t%M", "-o", str(metrics),
        str(binary), *arguments, "-o", str(output),
    ]
    started = time.perf_counter()
    with stdout.open("w", encoding="utf-8") as stdout_handle, stderr.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            command, stdout=stdout_handle, stderr=stderr_handle, check=False
        )
    wall_seconds = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {completed.returncode}: {' '.join(command)}"
        )
    elapsed, user, system, max_rss = metrics.read_text(
        encoding="utf-8"
    ).strip().split("\t")
    return {
        "wall_seconds": wall_seconds,
        "time_elapsed_seconds": float(elapsed),
        "user_seconds": float(user),
        "system_seconds": float(system),
        "max_rss_kib": int(max_rss),
        "network_data_sha256": network_data_hash(output),
    }


def median(rows, field):
    return statistics.median(row[field] for row in rows)


def main():
    args = parse_args()
    if args.repetitions < 1:
        raise SystemExit("--repetitions must be positive")

    for path in (args.baseline, args.candidate):
        if not path.is_file():
            raise SystemExit(f"Executable not found: {path}")

    args.output.mkdir(parents=True, exist_ok=True)
    empty_adjacency = args.output / "empty.adj"
    empty_adjacency.write_text("", encoding="utf-8")
    implementations = {"baseline": args.baseline, "candidate": args.candidate}
    benchmark_cases = cases(args.rank_data, args.brca_data, empty_adjacency)

    rows = []
    for case_name, arguments in benchmark_cases.items():
        for implementation, binary in implementations.items():
            prefix = args.output / f"warmup_{case_name}_{implementation}"
            run_once(
                binary, arguments, prefix.with_suffix(".adj"),
                prefix.with_suffix(".time"), prefix.with_suffix(".stdout"),
                prefix.with_suffix(".stderr"),
            )

        for repetition in range(1, args.repetitions + 1):
            order = (
                ("baseline", "candidate")
                if repetition % 2 else ("candidate", "baseline")
            )
            for implementation in order:
                prefix = args.output / (
                    f"{case_name}_{implementation}_r{repetition:02d}"
                )
                result = run_once(
                    implementations[implementation], arguments,
                    prefix.with_suffix(".adj"), prefix.with_suffix(".time"),
                    prefix.with_suffix(".stdout"), prefix.with_suffix(".stderr"),
                )
                rows.append(
                    {
                        "case": case_name,
                        "implementation": implementation,
                        "repetition": repetition,
                        **result,
                    }
                )

    raw_path = args.output / "raw_timings.csv"
    with raw_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    for case_name in benchmark_cases:
        by_implementation = {
            implementation: [
                row for row in rows
                if row["case"] == case_name
                and row["implementation"] == implementation
            ]
            for implementation in implementations
        }
        hashes = {
            row["network_data_sha256"]
            for implementation_rows in by_implementation.values()
            for row in implementation_rows
        }
        if len(hashes) != 1:
            raise RuntimeError(f"Network-data mismatch for {case_name}: {hashes}")

        baseline_wall = median(by_implementation["baseline"], "wall_seconds")
        candidate_wall = median(by_implementation["candidate"], "wall_seconds")
        baseline_rss = median(by_implementation["baseline"], "max_rss_kib")
        candidate_rss = median(by_implementation["candidate"], "max_rss_kib")
        summary_rows.append(
            {
                "case": case_name,
                "baseline_wall_seconds": baseline_wall,
                "candidate_wall_seconds": candidate_wall,
                "wall_speedup": baseline_wall / candidate_wall,
                "seconds_saved": baseline_wall - candidate_wall,
                "baseline_max_rss_kib": baseline_rss,
                "candidate_max_rss_kib": candidate_rss,
                "rss_ratio": candidate_rss / baseline_rss,
                "network_data_sha256": hashes.pop(),
            }
        )

    summary_path = args.output / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys())
        writer.writeheader()
        writer.writerows(summary_rows)

    print(summary_path.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
