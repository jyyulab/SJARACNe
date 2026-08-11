#!/usr/bin/env python3

"""Generate sparse multi-hub and dense DPI benchmark adjacency files."""

import argparse
import csv
from pathlib import Path


HUB_COUNTS = (10, 50, 100)
TARGET_COUNT = 4000
DENSE_GENE_COUNT = 500


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_accessions(path):
    with path.open(encoding="utf-8") as handle:
        next(handle)
        accessions = [line.split("\t", 1)[0] for line in handle if line.strip()]

    if len(accessions) < 5000:
        raise ValueError("This benchmark requires at least 5,000 genes")
    if len(accessions) != len(set(accessions)):
        raise ValueError("Expression accessions must be unique")
    return accessions


def target_mi(position):
    return f"{0.6 + 0.399 * (TARGET_COUNT - position) / TARGET_COUNT:.12g}"


def write_hub_case(output_dir, accessions, hub_count):
    hubs = accessions[:hub_count]
    targets = accessions[500 : 500 + TARGET_COUNT]
    adjacency = output_dir / f"hub_h{hub_count:04d}_t{TARGET_COUNT:04d}.adj"
    hub_file = output_dir / f"hub_h{hub_count:04d}.txt"

    hub_file.write_text("".join(f"{hub}\n" for hub in hubs), encoding="utf-8")
    with adjacency.open("w", encoding="utf-8", newline="\n") as handle:
        for hub_index, hub in enumerate(hubs):
            fields = [hub]
            for position, target in enumerate(targets):
                fields.extend((target, target_mi(position)))

            ring = {
                (hub_index - 2) % hub_count,
                (hub_index - 1) % hub_count,
                (hub_index + 1) % hub_count,
                (hub_index + 2) % hub_count,
            }
            for neighbor in sorted(ring):
                fields.extend((hubs[neighbor], "0.15"))
            handle.write("\t".join(fields) + "\n")

    return {
        "case_id": f"hub_h{hub_count:04d}_t{TARGET_COUNT:04d}",
        "mode": "selected_hubs",
        "adjacency": adjacency.name,
        "hubs": hub_file.name,
        "source_rows": hub_count,
        "genes_in_graph": hub_count + TARGET_COUNT,
        "directed_edges": hub_count * (TARGET_COUNT + 4),
    }


def write_all_gene_skew_case(output_dir, accessions):
    hub_count = 100
    hubs = accessions[:hub_count]
    targets = accessions[500 : 500 + TARGET_COUNT]
    adjacency = output_dir / "allgene_skew_h0100_t4000.adj"

    with adjacency.open("w", encoding="utf-8", newline="\n") as handle:
        for hub in hubs:
            fields = [hub]
            for position, target in enumerate(targets):
                fields.extend((target, target_mi(position)))
            handle.write("\t".join(fields) + "\n")

        for position, target in enumerate(targets):
            fields = [target]
            mi = target_mi(position)
            for hub in hubs:
                fields.extend((hub, mi))
            handle.write("\t".join(fields) + "\n")

    return {
        "case_id": "allgene_skew_h0100_t4000",
        "mode": "all_gene_skewed",
        "adjacency": adjacency.name,
        "hubs": "",
        "source_rows": hub_count + TARGET_COUNT,
        "genes_in_graph": hub_count + TARGET_COUNT,
        "directed_edges": 2 * hub_count * TARGET_COUNT,
    }


def write_dense_tied_case(output_dir, accessions):
    genes = accessions[:DENSE_GENE_COUNT]
    adjacency = output_dir / f"allgene_dense_tied_g{DENSE_GENE_COUNT:04d}.adj"

    with adjacency.open("w", encoding="utf-8", newline="\n") as handle:
        for source in genes:
            fields = [source]
            for target in genes:
                if target != source:
                    fields.extend((target, "0.5"))
            handle.write("\t".join(fields) + "\n")

    return {
        "case_id": f"allgene_dense_tied_g{DENSE_GENE_COUNT:04d}",
        "mode": "all_gene_dense_control",
        "adjacency": adjacency.name,
        "hubs": "",
        "source_rows": DENSE_GENE_COUNT,
        "genes_in_graph": DENSE_GENE_COUNT,
        "directed_edges": DENSE_GENE_COUNT * (DENSE_GENE_COUNT - 1),
    }


def main():
    args = parse_args()
    accessions = read_accessions(args.expression)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cases = [
        write_hub_case(args.output_dir, accessions, hub_count)
        for hub_count in HUB_COUNTS
    ]
    cases.append(write_all_gene_skew_case(args.output_dir, accessions))
    cases.append(write_dense_tied_case(args.output_dir, accessions))

    manifest = args.output_dir / "cases.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    for case in cases:
        print(
            f"{case['case_id']}: rows={case['source_rows']}, "
            f"edges={case['directed_edges']}"
        )


if __name__ == "__main__":
    main()
