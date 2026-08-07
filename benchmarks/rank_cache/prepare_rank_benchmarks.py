#!/usr/bin/env python3
"""Create deterministic expression panels for the AP-MI rank-cache benchmark."""

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
from contextlib import ExitStack
from pathlib import Path


HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[1]
DEFAULT_SPEC = HERE / "spec.json"
DEFAULT_TF_LIST = PROJECT_ROOT / "SJARACNe" / "config" / "TF_list.txt"
DEFAULT_OUTPUT = HERE / "generated"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values):
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def write_lf(path, text):
    with Path(path).open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def stable_order(values, salt, namespace):
    def key(value):
        payload = "{}\0{}\0{}".format(salt, namespace, value).encode("utf-8")
        return hashlib.sha256(payload).digest(), value

    return sorted(values, key=key)


def load_identifier_list(path):
    identifiers = []
    seen = set()
    duplicates = 0
    with Path(path).open("r", encoding="utf-8-sig", newline=None) as handle:
        for line in handle:
            identifier = line.strip()
            if not identifier:
                continue
            if identifier in seen:
                duplicates += 1
                continue
            seen.add(identifier)
            identifiers.append(identifier)
    return identifiers, duplicates


def inspect_expression(path, metadata_columns):
    genes = []
    seen_genes = set()
    with Path(path).open("r", encoding="utf-8-sig", newline=None) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError("Expression file is empty: {}".format(path))
        header = header_line.rstrip("\r\n").split("\t")
        if len(header) <= metadata_columns:
            raise ValueError(
                "Expression header must contain more than {} columns".format(
                    metadata_columns
                )
            )
        observations = header[metadata_columns:]
        if len(observations) != len(set(observations)):
            raise ValueError("Expression header contains duplicate observation IDs")

        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\r\n").split("\t", metadata_columns)
            if len(fields) <= metadata_columns:
                raise ValueError(
                    "Expression row {} has too few columns".format(line_number)
                )
            raw_gene = fields[0]
            gene = raw_gene.strip()
            if not gene:
                raise ValueError(
                    "Expression row {} has an empty first-column ID".format(
                        line_number
                    )
                )
            if raw_gene != gene:
                raise ValueError(
                    "Expression row {} has leading or trailing whitespace in its "
                    "first-column ID: {!r}".format(line_number, raw_gene)
                )
            if gene in seen_genes:
                raise ValueError("Duplicate expression gene ID: {}".format(gene))
            seen_genes.add(gene)
            genes.append(gene)

    return header[:metadata_columns], observations, genes


def resolve_count(value, total, label):
    if value == "all":
        return total
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("{} must be an integer or 'all'".format(label))
    if value < 1 or value > total:
        raise ValueError(
            "{}={} is outside the available range 1..{}".format(label, value, total)
        )
    return value


def required_unique_strings(spec, key):
    values = spec.get(key, [])
    if not isinstance(values, list) or any(
        not isinstance(value, str) or not value for value in values
    ):
        raise ValueError("{} must be a list of non-empty strings".format(key))
    if len(values) != len(set(values)):
        raise ValueError("{} contains duplicates".format(key))
    return values


def validate_expected_source(spec, actual):
    expected = spec.get("expected_source")
    if expected is None:
        return
    if not isinstance(expected, dict):
        raise ValueError("expected_source must be an object")
    for key, expected_value in expected.items():
        if key not in actual:
            raise ValueError("Unknown expected_source field: {}".format(key))
        if actual[key] != expected_value:
            raise ValueError(
                "Source identity mismatch for {}: expected {!r}, found {!r}".format(
                    key, expected_value, actual[key]
                )
            )


def build_gene_and_hub_orders(genes, tf_ids, spec):
    gene_set = set(genes)
    eligible_hubs = gene_set.intersection(tf_ids)
    required_hubs = required_unique_strings(spec, "required_hubs")
    required_targets = required_unique_strings(spec, "required_targets")

    missing_hubs = [gene for gene in required_hubs if gene not in eligible_hubs]
    if missing_hubs:
        raise ValueError(
            "Required hub genes are absent from the expression/TF overlap: {}".format(
                ", ".join(missing_hubs)
            )
        )
    missing_targets = [gene for gene in required_targets if gene not in gene_set]
    if missing_targets:
        raise ValueError(
            "Required target genes are absent from the expression matrix: {}".format(
                ", ".join(missing_targets)
            )
        )
    target_hub_overlap = set(required_targets).intersection(eligible_hubs)
    if target_hub_overlap:
        raise ValueError(
            "Required targets cannot also be eligible hubs: {}".format(
                ", ".join(sorted(target_hub_overlap))
            )
        )

    salt = spec["selection_salt"]
    hub_order = required_hubs + stable_order(
        eligible_hubs.difference(required_hubs), salt, "hub"
    )
    non_hubs = gene_set.difference(eligible_hubs)
    non_hub_order = required_targets + stable_order(
        non_hubs.difference(required_targets), salt, "target"
    )

    raw_gene_counts = spec["gene_sweep"]["genes"]
    numeric_gene_counts = [value for value in raw_gene_counts if value != "all"]
    if not numeric_gene_counts:
        raise ValueError("gene_sweep.genes must contain at least one integer")
    smallest_gene_panel = min(numeric_gene_counts)
    core_hubs = spec["core_hubs_in_smallest_gene_panel"]
    if isinstance(core_hubs, bool) or not isinstance(core_hubs, int):
        raise ValueError("core_hubs_in_smallest_gene_panel must be an integer")
    if core_hubs < 1 or core_hubs > len(hub_order):
        raise ValueError(
            "core_hubs_in_smallest_gene_panel={} is outside 1..{}".format(
                core_hubs, len(hub_order)
            )
        )
    non_hubs_in_prefix = smallest_gene_panel - core_hubs
    if non_hubs_in_prefix < len(required_targets):
        raise ValueError(
            "The smallest gene panel cannot contain the core hubs and required targets"
        )
    if non_hubs_in_prefix > len(non_hub_order):
        raise ValueError("Not enough non-hub genes to fill the smallest gene panel")

    master_gene_order = (
        hub_order[:core_hubs]
        + non_hub_order[:non_hubs_in_prefix]
        + hub_order[core_hubs:]
        + non_hub_order[non_hubs_in_prefix:]
    )
    if len(master_gene_order) != len(genes) or len(set(master_gene_order)) != len(
        genes
    ):
        raise AssertionError("Internal error: master gene order is not a permutation")
    return hub_order, master_gene_order


def build_cases(spec, gene_total, observation_total, hub_total):
    cases = {}

    def add_case(sweep, gene_value, observation_value, hub_value):
        genes = resolve_count(gene_value, gene_total, "genes")
        observations = resolve_count(
            observation_value, observation_total, "observations"
        )
        hubs = resolve_count(hub_value, hub_total, "hubs")
        key = (genes, observations, hubs)
        cases.setdefault(key, set()).add(sweep)

    hub_sweep = spec["hub_sweep"]
    for hubs in hub_sweep["hubs"]:
        add_case("hub_count", hub_sweep["genes"], hub_sweep["observations"], hubs)

    observation_sweep = spec["observation_sweep"]
    for observations in observation_sweep["observations"]:
        add_case(
            "observation_count",
            observation_sweep["genes"],
            observations,
            observation_sweep["hubs"],
        )

    gene_sweep = spec["gene_sweep"]
    for genes in gene_sweep["genes"]:
        add_case(
            "gene_count", genes, gene_sweep["observations"], gene_sweep["hubs"]
        )
    return cases


def expression_name(genes, observations):
    return "expression_g{:05d}_n{:04d}.exp".format(genes, observations)


def hub_name(hubs):
    return "hubs_h{:04d}.txt".format(hubs)


def case_name(genes, observations, hubs):
    return "g{:05d}_n{:04d}_h{:04d}".format(genes, observations, hubs)


def write_expression_panels(
    source,
    staging,
    metadata_columns,
    metadata_header,
    observations,
    observation_order,
    master_gene_order,
    matrix_keys,
):
    gene_panels = {genes: set(master_gene_order[:genes]) for genes, _ in matrix_keys}
    cell_indices = {
        count: observation_order[:count]
        for count in {observations_count for _, observations_count in matrix_keys}
    }
    written = {key: 0 for key in matrix_keys}

    with ExitStack() as stack:
        outputs = {}
        for key in matrix_keys:
            genes, observation_count = key
            output_path = staging / expression_name(genes, observation_count)
            handle = stack.enter_context(
                output_path.open("w", encoding="utf-8", newline="\n")
            )
            selected_indices = cell_indices[observation_count]
            selected_headers = [observations[index] for index in selected_indices]
            handle.write("\t".join(metadata_header + selected_headers) + "\n")
            outputs[key] = (handle, selected_indices)

        with Path(source).open("r", encoding="utf-8-sig", newline=None) as handle:
            handle.readline()
            expected_columns = metadata_columns + len(observations)
            for line_number, line in enumerate(handle, start=2):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != expected_columns:
                    raise ValueError(
                        "Expression row {} has {} columns; expected {}".format(
                            line_number, len(fields), expected_columns
                        )
                    )
                for column_number, value in enumerate(
                    fields[metadata_columns:], start=metadata_columns + 1
                ):
                    try:
                        numeric = float(value)
                    except ValueError:
                        raise ValueError(
                            "Expression row {}, column {} is not numeric: {!r}".format(
                                line_number, column_number, value
                            )
                        )
                    if not math.isfinite(numeric):
                        raise ValueError(
                            "Expression row {}, column {} is non-finite: {!r}".format(
                                line_number, column_number, value
                            )
                        )
                gene = fields[0].strip()
                for key, (output, selected_indices) in outputs.items():
                    genes, _ = key
                    if gene not in gene_panels[genes]:
                        continue
                    values = [fields[metadata_columns + index] for index in selected_indices]
                    output.write("\t".join(fields[:metadata_columns] + values) + "\n")
                    written[key] += 1

    for key, count in written.items():
        if count != key[0]:
            raise AssertionError(
                "Internal error: wrote {} genes for {}".format(count, key)
            )


def artifact_record(path, kind, genes=None, observations=None, hubs=None):
    record = {
        "path": path.name,
        "kind": kind,
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if genes is not None:
        record["genes"] = genes
    if observations is not None:
        record["observations"] = observations
    if hubs is not None:
        record["hubs"] = hubs
    return record


def prepare(expression, tf_list, spec_path, output_dir, force=False):
    expression = Path(expression).resolve()
    tf_list = Path(tf_list).resolve()
    spec_path = Path(spec_path).resolve()
    output_dir = Path(output_dir).resolve()
    if not expression.is_file():
        raise ValueError("Expression file does not exist: {}".format(expression))
    if not tf_list.is_file():
        raise ValueError("TF list does not exist: {}".format(tf_list))

    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    if spec.get("schema_version") != 1:
        raise ValueError("Unsupported benchmark specification schema")
    if not isinstance(spec.get("selection_salt"), str) or not spec["selection_salt"]:
        raise ValueError("selection_salt must be a non-empty string")
    metadata_columns = spec.get("metadata_columns")
    if isinstance(metadata_columns, bool) or not isinstance(metadata_columns, int):
        raise ValueError("metadata_columns must be an integer")
    if metadata_columns < 1:
        raise ValueError("metadata_columns must be positive")

    print("Inspecting expression matrix: {}".format(expression))
    metadata_header, observations, genes = inspect_expression(
        expression, metadata_columns
    )
    expected_metadata_header = spec.get("expected_metadata_header")
    if (
        expected_metadata_header is not None
        and metadata_header != expected_metadata_header
    ):
        raise ValueError(
            "Metadata header mismatch: expected {!r}, found {!r}".format(
                expected_metadata_header, metadata_header
            )
        )
    tf_ids, tf_duplicates = load_identifier_list(tf_list)
    hub_order, master_gene_order = build_gene_and_hub_orders(genes, tf_ids, spec)
    expression_sha256 = sha256_file(expression)
    tf_list_sha256 = sha256_file(tf_list)
    tf_ids_sha256 = sha256_lines(tf_ids)
    validate_expected_source(
        spec,
        {
            "expression_sha256": expression_sha256,
            "tf_ids_sha256": tf_ids_sha256,
            "gene_count": len(genes),
            "observation_count": len(observations),
            "eligible_hubs": len(hub_order),
        },
    )
    ordered_observations = stable_order(
        observations, spec["selection_salt"], "observation"
    )
    observation_index = {value: index for index, value in enumerate(observations)}
    observation_order = [observation_index[value] for value in ordered_observations]
    cases = build_cases(spec, len(genes), len(observations), len(hub_order))

    matrix_keys = sorted(
        {(genes_count, obs_count) for genes_count, obs_count, _ in cases}
    )
    hub_counts = sorted({hubs for _, _, hubs in cases})
    for genes_count, _, hubs in cases:
        selected_genes = set(master_gene_order[:genes_count])
        missing = set(hub_order[:hubs]).difference(selected_genes)
        if missing:
            raise ValueError(
                "The g={} panel is missing {} requested hubs; increase the fixed gene "
                "count or revise the panel layout".format(genes_count, len(missing))
            )

    expected_names = {
        expression_name(genes_count, obs_count)
        for genes_count, obs_count in matrix_keys
    }
    expected_names.update(hub_name(hubs) for hubs in hub_counts)
    expected_names.update({"benchmark_cases.csv", "manifest.json"})
    existing = [
        output_dir / name
        for name in sorted(expected_names)
        if (output_dir / name).exists()
    ]
    if existing and not force:
        raise ValueError(
            "Benchmark artifacts already exist (use --force to replace them): {}".format(
                ", ".join(str(path) for path in existing[:5])
            )
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".rank_cache_staging_", dir=str(output_dir.parent)
    ) as temporary:
        staging = Path(temporary)
        print(
            "Writing {} unique expression matrices for {} benchmark cases".format(
                len(matrix_keys), len(cases)
            )
        )
        write_expression_panels(
            expression,
            staging,
            metadata_columns,
            metadata_header,
            observations,
            observation_order,
            master_gene_order,
            matrix_keys,
        )

        for hubs in hub_counts:
            write_lf(
                staging / hub_name(hubs),
                "".join("{}\n".format(gene) for gene in hub_order[:hubs]),
            )

        case_rows = []
        for (genes_count, obs_count, hubs), sweeps in sorted(cases.items()):
            case_rows.append(
                {
                    "case_id": case_name(genes_count, obs_count, hubs),
                    "sweeps": ";".join(sorted(sweeps)),
                    "expression": expression_name(genes_count, obs_count),
                    "hubs": hub_name(hubs),
                    "gene_count": genes_count,
                    "observation_count": obs_count,
                    "hub_count": hubs,
                    "candidate_mi_pairs": hubs * (genes_count - 1),
                }
            )
        with (staging / "benchmark_cases.csv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(case_rows[0]))
            writer.writeheader()
            writer.writerows(case_rows)

        artifacts = []
        for genes_count, obs_count in matrix_keys:
            artifacts.append(
                artifact_record(
                    staging / expression_name(genes_count, obs_count),
                    "expression",
                    genes=genes_count,
                    observations=obs_count,
                )
            )
        for hubs in hub_counts:
            artifacts.append(
                artifact_record(staging / hub_name(hubs), "hub_list", hubs=hubs)
            )
        artifacts.append(
            artifact_record(staging / "benchmark_cases.csv", "case_table")
        )

        manifest = {
            "schema_version": 1,
            "purpose": "adaptive-partitioning MI rank-cache performance benchmark",
            "source": {
                "expression_path": str(expression),
                "expression_bytes": expression.stat().st_size,
                "expression_sha256": expression_sha256,
                "gene_count": len(genes),
                "observation_count": len(observations),
                "tf_list_path": str(tf_list),
                "tf_list_sha256": tf_list_sha256,
                "tf_ids_sha256": tf_ids_sha256,
                "tf_ids": len(tf_ids),
                "tf_duplicates_ignored": tf_duplicates,
                "eligible_hubs": len(hub_order),
            },
            "specification": {
                "path": str(spec_path),
                "sha256": sha256_file(spec_path),
                "selection_salt": spec["selection_salt"],
            },
            "selection": {
                "gene_panels_are_nested": True,
                "observation_panels_are_nested": True,
                "output_gene_rows_follow_source_order": True,
                "hub_order_sha256": sha256_lines(hub_order),
                "master_gene_order_sha256": sha256_lines(master_gene_order),
                "observation_order_sha256": sha256_lines(ordered_observations),
                "required_hubs": spec.get("required_hubs", []),
                "required_targets": spec.get("required_targets", []),
            },
            "case_count": len(case_rows),
            "artifacts": artifacts,
        }
        write_lf(
            staging / "manifest.json",
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        )

        output_dir.mkdir(parents=True, exist_ok=True)
        for name in sorted(expected_names):
            os.replace(str(staging / name), str(output_dir / name))

    print("Created benchmark artifacts in {}".format(output_dir))
    print(
        "Source: {} genes x {} observations; eligible hubs: {}".format(
            len(genes), len(observations), len(hub_order)
        )
    )
    return output_dir / "manifest.json"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True, help="Source expression TSV")
    parser.add_argument(
        "--tf-list", default=str(DEFAULT_TF_LIST), help="Eligible TF identifiers"
    )
    parser.add_argument(
        "--spec", default=str(DEFAULT_SPEC), help="Benchmark JSON specification"
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT),
        help="Directory for generated benchmark artifacts",
    )
    parser.add_argument(
        "--force", action="store_true", help="Replace known generated artifacts"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        prepare(
            args.expression,
            args.tf_list,
            args.spec,
            args.output_dir,
            force=args.force,
        )
    except (OSError, ValueError, KeyError, TypeError) as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
