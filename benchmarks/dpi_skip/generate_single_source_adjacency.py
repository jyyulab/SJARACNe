#!/usr/bin/env python3

"""Generate a one-source adjacency row that exercises failed DPI lookups."""

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--hub-output", type=Path)
    args = parser.parse_args()

    with args.expression.open(encoding="utf-8") as handle:
        next(handle)
        accessions = [line.split("\t", 1)[0] for line in handle if line.strip()]

    if len(accessions) < 3:
        raise ValueError("At least three genes are required")
    if len(accessions) != len(set(accessions)):
        raise ValueError("Expression accessions must be unique")

    # Strictly decreasing positive MI values force the legacy DPI loop to test
    # every preceding neighbor pair instead of taking its equal-MI early exit.
    denominator = len(accessions)
    fields = [accessions[0]]
    for position, target in enumerate(accessions[1:], start=1):
        mi = (denominator - position) / denominator
        fields.extend((target, f"{mi:.12g}"))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\t".join(fields) + "\n", encoding="utf-8")
    if args.hub_output is not None:
        args.hub_output.parent.mkdir(parents=True, exist_ok=True)
        args.hub_output.write_text(accessions[0] + "\n", encoding="utf-8")
    print(
        f"wrote {args.output}: source_rows=1, "
        f"neighbors={len(accessions) - 1}"
    )


if __name__ == "__main__":
    main()
