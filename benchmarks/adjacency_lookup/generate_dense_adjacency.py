#!/usr/bin/env python3

"""Generate a deterministic dense adjacency fixture from an expression matrix."""

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sources", required=True, type=int)
    parser.add_argument("--mi", default="0.5")
    return parser.parse_args()


def read_accessions(path):
    with path.open(encoding="utf-8") as handle:
        next(handle)
        accessions = [line.split("\t", 1)[0] for line in handle if line.strip()]

    if len(accessions) != len(set(accessions)):
        raise ValueError("Expression accessions must be unique for this benchmark")
    return accessions


def main():
    args = parse_args()
    accessions = read_accessions(args.expression)

    if not 1 <= args.sources <= len(accessions):
        raise ValueError("--sources must be between 1 and the number of genes")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        for source in accessions[: args.sources]:
            fields = [source]
            for target in accessions:
                if target != source:
                    fields.extend((target, args.mi))
            handle.write("\t".join(fields) + "\n")

    edge_count = args.sources * (len(accessions) - 1)
    print(
        f"wrote {args.output}: genes={len(accessions)}, "
        f"sources={args.sources}, edges={edge_count}"
    )


if __name__ == "__main__":
    main()
