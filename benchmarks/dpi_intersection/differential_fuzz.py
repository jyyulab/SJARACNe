#!/usr/bin/env python3

"""Compare legacy and candidate DPI outputs on deterministic random networks."""

import argparse
import random
import subprocess
import tempfile
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--cases", type=int, default=250)
    parser.add_argument("--seed", type=int, default=20260808)
    return parser.parse_args()


def data_rows(path):
    return b"".join(
        line
        for line in path.read_bytes().splitlines(keepends=True)
        if line.strip() and not line.startswith(b">")
    )


def write_expression(path, genes):
    lines = ["isoformId\tgeneSymbol\ts1\ts2\ts3\ts4"]
    for index, gene in enumerate(genes):
        values = [str((index * multiplier) % 11 + 1) for multiplier in (1, 3, 5, 7)]
        lines.append("\t".join((gene, gene, *values)))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_adjacency(path, genes, sources, rng, symmetric):
    mi_values = ("0.1", "0.2", "0.3", "0.5", "0.7", "0.9")
    edges = {}

    if symmetric:
        for left in range(len(genes)):
            for right in range(left + 1, len(genes)):
                if rng.random() < 0.35:
                    mi = rng.choice(mi_values)
                    edges[left, right] = mi
                    edges[right, left] = mi
    else:
        for source in sources:
            for target in range(len(genes)):
                if source != target and rng.random() < 0.35:
                    edges[source, target] = rng.choice(mi_values)

    lines = []
    for source in sources:
        pairs = [
            (target, edges[source, target])
            for target in range(len(genes))
            if (source, target) in edges
        ]

        # Exercise the adjacency representation rather than feeding it only
        # canonical rows: targets can be unordered, a target can occur twice,
        # and one source can be split across repeated lines. The final retained
        # occurrence remains the defined value for both implementations.
        if pairs and rng.random() < 0.4:
            duplicate_target, _ = rng.choice(pairs)
            pairs.append((duplicate_target, rng.choice(mi_values)))
        rng.shuffle(pairs)

        if not pairs:
            lines.append(genes[source])
            continue

        record_count = rng.randint(1, min(3, len(pairs)))
        records = [[] for _ in range(record_count)]
        for index, pair in enumerate(pairs):
            records[index % record_count].append(pair)

        for record in records:
            fields = [genes[source]]
            for target, mi in record:
                fields.extend((genes[target], mi))
            lines.append("\t".join(fields))

    rng.shuffle(lines)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(executable, expression, adjacency, output, epsilon, hubs, tf_genes):
    command = [
        str(executable),
        "-i",
        str(expression),
        "-j",
        str(adjacency),
        "-p",
        "1",
        "-e",
        str(epsilon),
        "-o",
        str(output),
    ]
    if hubs is not None:
        command.extend(("-s", str(hubs)))
    if tf_genes is not None:
        command.extend(("-l", str(tf_genes)))

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(result.stdout + result.stderr)


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    with tempfile.TemporaryDirectory() as folder_name:
        folder = Path(folder_name)

        for case in range(args.cases):
            gene_count = rng.randint(6, 16)
            genes = [f"G{index:02d}" for index in range(gene_count)]
            symmetric = rng.random() < 0.4

            if symmetric:
                sources = list(range(gene_count))
            else:
                source_count = rng.randint(2, max(2, gene_count // 2))
                sources = sorted(rng.sample(range(gene_count), source_count))

            expression = folder / "input.exp"
            adjacency = folder / "input.adj"
            write_expression(expression, genes)
            write_adjacency(adjacency, genes, sources, rng, symmetric)

            hubs = None
            if not symmetric and rng.random() < 0.8:
                selected = rng.sample(sources, rng.randint(1, len(sources)))
                hubs = folder / "hubs.txt"
                hubs.write_text(
                    "".join(f"{genes[index]}\n" for index in selected),
                    encoding="utf-8",
                )

            tf_genes = None
            if rng.random() < 0.5:
                selected = [gene for gene in genes if rng.random() < 0.3]
                if selected:
                    tf_genes = folder / "tf_genes.txt"
                    tf_genes.write_text(
                        "".join(f"{gene}\n" for gene in selected),
                        encoding="utf-8",
                    )

            epsilon = rng.choice((0, 0.1, 0.25, 0.5, 0.9))
            baseline_output = folder / "baseline.adj"
            candidate_output = folder / "candidate.adj"
            run(
                args.baseline,
                expression,
                adjacency,
                baseline_output,
                epsilon,
                hubs,
                tf_genes,
            )
            run(
                args.candidate,
                expression,
                adjacency,
                candidate_output,
                epsilon,
                hubs,
                tf_genes,
            )

            if data_rows(baseline_output) != data_rows(candidate_output):
                raise RuntimeError(
                    f"DPI output mismatch in case {case} (seed={args.seed})"
                )

    print(
        f"matched {args.cases} differential cases "
        f"(seed={args.seed})"
    )


if __name__ == "__main__":
    main()
