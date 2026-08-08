#!/usr/bin/env python3

"""Run matched legacy/candidate DPI intersection benchmarks."""

import argparse
import csv
import hashlib
import statistics
import subprocess
import time
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--fixtures-dir", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    return parser.parse_args()


def data_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for line in handle:
            if line.strip() and not line.startswith(b">"):
                digest.update(line)
    return digest.hexdigest()


def main():
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    with (args.fixtures_dir / "cases.csv").open(encoding="utf-8", newline="") as handle:
        cases = list(csv.DictReader(handle))

    implementations = {
        "baseline": args.baseline,
        "candidate": args.candidate,
    }
    records = []

    for case in cases:
        adjacency = args.fixtures_dir / case["adjacency"]
        hubs = args.fixtures_dir / case["hubs"] if case["hubs"] else None

        for repeat in range(1, args.repeats + 1):
            order = ("baseline", "candidate")
            if repeat % 2 == 0:
                order = tuple(reversed(order))

            for implementation in order:
                output = args.results_dir / (
                    f"{case['case_id']}_{implementation}_r{repeat}.adj"
                )
                command = [
                    str(implementations[implementation]),
                    "-i",
                    str(args.expression),
                    "-j",
                    str(adjacency),
                    "-p",
                    "1",
                    "-e",
                    "0",
                    "-o",
                    str(output),
                ]
                if hubs is not None:
                    command.extend(("-s", str(hubs)))

                started = time.perf_counter()
                result = subprocess.run(command, capture_output=True, text=True)
                elapsed = time.perf_counter() - started

                if result.returncode != 0:
                    raise RuntimeError(result.stdout + result.stderr)
                if "[NETWORK] Applying DPI" not in result.stdout:
                    raise RuntimeError(
                        f"DPI did not run for {case['case_id']} {implementation}"
                    )

                digest = data_sha256(output)
                records.append(
                    {
                        **case,
                        "implementation": implementation,
                        "repeat": repeat,
                        "wall_seconds": f"{elapsed:.6f}",
                        "data_sha256": digest,
                    }
                )
                print(
                    f"{case['case_id']} {implementation} r{repeat}: "
                    f"{elapsed:.3f} s",
                    flush=True,
                )
                output.unlink()

    for case in cases:
        hashes = {
            row["data_sha256"]
            for row in records
            if row["case_id"] == case["case_id"]
        }
        if len(hashes) != 1:
            raise RuntimeError(
                f"Network outputs differ for {case['case_id']}: {sorted(hashes)}"
            )

    raw_path = args.results_dir / "raw_timings.csv"
    with raw_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    summaries = []
    for case in cases:
        times = {}
        for implementation in implementations:
            times[implementation] = [
                float(row["wall_seconds"])
                for row in records
                if row["case_id"] == case["case_id"]
                and row["implementation"] == implementation
            ]
        baseline = statistics.median(times["baseline"])
        candidate = statistics.median(times["candidate"])
        summaries.append(
            {
                **case,
                "baseline_median_seconds": f"{baseline:.6f}",
                "candidate_median_seconds": f"{candidate:.6f}",
                "speedup": f"{baseline / candidate:.6f}",
            }
        )

    with (args.results_dir / "summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)


if __name__ == "__main__":
    main()
