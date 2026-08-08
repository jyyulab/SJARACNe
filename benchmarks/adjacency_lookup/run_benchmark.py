#!/usr/bin/env python3

"""Benchmark native -j adjacency loading with matched baseline/candidate builds."""

import argparse
import csv
import hashlib
import subprocess
import time
from pathlib import Path


CASES = (
    ("g01000_s0100", 1000, 100, 99900),
    ("g05000_s0020", 5000, 20, 99980),
    ("g19936_s0005", 19936, 5, 99675),
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--expression-dir", required=True, type=Path)
    parser.add_argument("--adjacency-dir", required=True, type=Path)
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
    records = []

    implementations = {
        "baseline": args.baseline,
        "candidate": args.candidate,
    }

    for case_id, genes, sources, edges in CASES:
        expression = args.expression_dir / f"expression_g{genes:05d}_n0500.exp"
        adjacency = args.adjacency_dir / f"{case_id}.adj"

        for repeat in range(1, args.repeats + 1):
            order = ("baseline", "candidate")
            if repeat % 2 == 0:
                order = tuple(reversed(order))

            for implementation in order:
                output = args.results_dir / f"{case_id}_{implementation}_r{repeat}.adj"
                command = [
                    str(implementations[implementation]),
                    "-i",
                    str(expression),
                    "-j",
                    str(adjacency),
                    "-p",
                    "1",
                    "-e",
                    "1",
                    "-o",
                    str(output),
                ]

                started = time.perf_counter()
                result = subprocess.run(command, capture_output=True, text=True)
                elapsed = time.perf_counter() - started

                stdout_path = output.with_suffix(".stdout.txt")
                stderr_path = output.with_suffix(".stderr.txt")
                stdout_path.write_text(result.stdout, encoding="utf-8")
                stderr_path.write_text(result.stderr, encoding="utf-8")

                if result.returncode != 0:
                    raise RuntimeError(
                        f"{case_id} {implementation} repeat {repeat} failed: "
                        f"{result.stderr}"
                    )

                record = {
                    "case_id": case_id,
                    "genes": genes,
                    "sources": sources,
                    "retained_edges": edges,
                    "implementation": implementation,
                    "repeat": repeat,
                    "wall_seconds": f"{elapsed:.6f}",
                    "data_sha256": data_sha256(output),
                }
                records.append(record)
                print(
                    f"{case_id} {implementation} r{repeat}: {elapsed:.3f} s",
                    flush=True,
                )

    for case_id, _, _, _ in CASES:
        hashes = {
            record["data_sha256"]
            for record in records
            if record["case_id"] == case_id
        }
        if len(hashes) != 1:
            raise RuntimeError(f"Output data differ for {case_id}: {sorted(hashes)}")

    fieldnames = list(records[0])
    with (args.results_dir / "raw_timings.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
