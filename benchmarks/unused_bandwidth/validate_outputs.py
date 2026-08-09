#!/usr/bin/env python3

import argparse
import hashlib
import subprocess
import tempfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare retained network data before and after bandwidth removal."
    )
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--brca-data", type=Path, required=True)
    return parser.parse_args()


def network_data(path):
    return b"".join(
        line
        for line in path.read_bytes().splitlines(keepends=True)
        if line.strip() and not line.startswith(b">")
    )


def run(binary, arguments, output):
    completed = subprocess.run(
        [str(binary), *arguments, "-o", str(output)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{binary} failed with exit code {completed.returncode}:\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    return network_data(output)


def main():
    args = parse_args()
    expression = args.fixtures / "tied_counts.exp"
    hubs = args.fixtures / "tied_hubs.txt"
    common = ["-i", str(expression), "-t", "0", "-e", "1"]
    cases = {
        "bootstrap": [*common, "-s", str(hubs), "-S", "17", "-r", "1"],
        "conditional_bootstrap": [
            *common, "-s", str(hubs), "-S", "17",
            "-c", "-_FOXP1", "0.75", "-r", "1",
        ],
        "noise_corrected_bootstrap": [
            *common, "-s", str(hubs), "-S", "17", "-r", "1", "-n", "0.1",
        ],
        "all_gene": [*common, "-S", "17", "-r", "0"],
        "brca100_adjacency_replay": [
            "-i", str(args.brca_data / "BRCA100.exp"),
            "-j", str(args.brca_data / "TF_run.adj"),
            "-s", str(args.brca_data / "BRCA100_TF.txt"),
            "-l", str(args.brca_data / "BRCA100_TF.txt"),
            "-p", "1", "-e", "0",
        ],
    }

    with tempfile.TemporaryDirectory(prefix="sjaracne-bandwidth-validation-") as folder:
        workdir = Path(folder)
        for name, arguments in cases.items():
            baseline_data = run(
                args.baseline, arguments, workdir / f"{name}_baseline.adj"
            )
            candidate_data = run(
                args.candidate, arguments, workdir / f"{name}_candidate.adj"
            )
            if not baseline_data:
                raise RuntimeError(f"Validation case produced no retained data: {name}")
            if candidate_data != baseline_data:
                raise RuntimeError(f"Network-data mismatch: {name}")
            digest = hashlib.sha256(baseline_data).hexdigest()
            print(f"{name}\t{digest}")


if __name__ == "__main__":
    main()
