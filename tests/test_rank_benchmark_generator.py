#!/usr/bin/env python3

import csv
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = (
    PROJECT_ROOT / "benchmarks" / "rank_cache" / "prepare_rank_benchmarks.py"
)
FIXTURES = PROJECT_ROOT / "benchmarks" / "rank_cache" / "fixtures"


def find_sjaracne_executable():
    configured = os.environ.get("SJARACNE_TEST_EXE")
    candidates = [
        configured,
        PROJECT_ROOT / "SJARACNe" / "bin" / "sjaracne.exe",
        shutil.which("sjaracne.exe"),
    ]
    for candidate in candidates:
        if not candidate or not Path(candidate).is_file():
            continue
        if os.name == "nt":
            with Path(candidate).open("rb") as handle:
                if handle.read(2) != b"MZ":
                    continue
        return str(candidate)
    return None


SJARACNE_EXE = find_sjaracne_executable()
REFERENCE_RUNTIME_SUPPORTED = (
    sys.platform.startswith("linux") and platform.libc_ver()[0] == "glibc"
)


class TestRankBenchmarkGenerator(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.workdir = Path(self.temporary.name)
        self.expression = self.workdir / "source.exp"
        self.tf_list = self.workdir / "tf.txt"
        self.spec = self.workdir / "spec.json"
        self.output = self.workdir / "generated"

        cells = ["c{:02d}".format(index) for index in range(1, 13)]
        genes = ["FOXP1", "CDKN1A"] + [
            "G{:02d}".format(index) for index in range(1, 13)
        ]
        rows = ["isoformId\tgeneSymbol\t{}".format("\t".join(cells))]
        for gene_index, gene in enumerate(genes):
            values = [
                str((gene_index * 3 + cell_index) % 7)
                for cell_index in range(12)
            ]
            rows.append("{}\t{}\t{}".format(gene, gene, "\t".join(values)))
        self.expression.write_text("\n".join(rows) + "\n", encoding="utf-8")
        self.tf_list.write_bytes(
            b"\xef\xbb\xbfFOXP1\r\nG01\r\nG02\r\nG03\r\nG04\r\nG05\r\nG05\r\n"
        )
        self.spec.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "selection_salt": "unit-test",
                    "metadata_columns": 2,
                    "required_hubs": ["FOXP1"],
                    "required_targets": ["CDKN1A"],
                    "core_hubs_in_smallest_gene_panel": 2,
                    "hub_sweep": {
                        "genes": 10,
                        "observations": 12,
                        "hubs": [1, 2, 4, "all"],
                    },
                    "observation_sweep": {
                        "genes": 10,
                        "hubs": 2,
                        "observations": [4, 8, 12],
                    },
                    "gene_sweep": {
                        "observations": 8,
                        "hubs": 2,
                        "genes": [6, 10, "all"],
                    },
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.temporary.cleanup()

    def run_generator(self, *extra):
        return subprocess.run(
            [
                sys.executable,
                str(GENERATOR),
                "--expression",
                str(self.expression),
                "--tf-list",
                str(self.tf_list),
                "--spec",
                str(self.spec),
                "--output-dir",
                str(self.output),
                *extra,
            ],
            capture_output=True,
            text=True,
            check=False,
        )

    @staticmethod
    def expression_axes(path):
        with Path(path).open("r", encoding="utf-8") as handle:
            rows = [line.rstrip("\n").split("\t") for line in handle]
        return rows[0][2:], {row[0] for row in rows[1:]}

    def test_generation_is_nested_deterministic_and_fail_closed(self):
        first = self.run_generator()
        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)

        manifest = json.loads(
            (self.output / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["source"]["gene_count"], 14)
        self.assertEqual(manifest["source"]["observation_count"], 12)
        self.assertEqual(manifest["source"]["eligible_hubs"], 6)
        self.assertEqual(manifest["source"]["tf_duplicates_ignored"], 1)
        self.assertEqual(manifest["case_count"], 8)

        with (self.output / "benchmark_cases.csv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            cases = list(csv.DictReader(handle))
        self.assertEqual(len(cases), 8)
        largest_hub_case = [
            row for row in cases if row["case_id"] == "g00010_n0012_h0006"
        ]
        self.assertEqual(
            {int(row["candidate_mi_pairs"]) for row in largest_hub_case}, {54}
        )

        hub_files = [
            self.output / "hubs_h0001.txt",
            self.output / "hubs_h0002.txt",
            self.output / "hubs_h0004.txt",
            self.output / "hubs_h0006.txt",
        ]
        hub_lists = [
            path.read_text(encoding="utf-8").splitlines() for path in hub_files
        ]
        self.assertEqual(hub_lists[0], ["FOXP1"])
        for smaller, larger in zip(hub_lists, hub_lists[1:]):
            self.assertEqual(larger[: len(smaller)], smaller)

        cells4, _ = self.expression_axes(
            self.output / "expression_g00010_n0004.exp"
        )
        cells8, genes10 = self.expression_axes(
            self.output / "expression_g00010_n0008.exp"
        )
        cells12, _ = self.expression_axes(
            self.output / "expression_g00010_n0012.exp"
        )
        _, genes6 = self.expression_axes(
            self.output / "expression_g00006_n0008.exp"
        )
        _, genes14 = self.expression_axes(
            self.output / "expression_g00014_n0008.exp"
        )
        self.assertEqual(cells8[:4], cells4)
        self.assertEqual(cells12[:8], cells8)
        self.assertLess(genes6, genes10)
        self.assertLess(genes10, genes14)
        self.assertIn("CDKN1A", genes6)
        self.assertTrue(set(hub_lists[-1]).issubset(genes10))

        hashes_before = {
            item["path"]: item["sha256"] for item in manifest["artifacts"]
        }
        refused = self.run_generator()
        self.assertEqual(refused.returncode, 1)
        self.assertIn("already exist", refused.stderr)

        replaced = self.run_generator("--force")
        self.assertEqual(replaced.returncode, 0, replaced.stdout + replaced.stderr)
        regenerated = json.loads(
            (self.output / "manifest.json").read_text(encoding="utf-8")
        )
        hashes_after = {
            item["path"]: item["sha256"] for item in regenerated["artifacts"]
        }
        self.assertEqual(hashes_after, hashes_before)


@unittest.skipUnless(
    SJARACNE_EXE and REFERENCE_RUNTIME_SUPPORTED,
    "fixed bootstrap references require a Linux glibc build",
)
class TestRankBenchmarkCorrectnessFixture(unittest.TestCase):
    def assert_fixture_matches_reference(self, reference_name, *extra_args):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "network.adj"
            result = subprocess.run(
                [
                    SJARACNE_EXE,
                    "-i",
                    str(FIXTURES / "tied_counts.exp"),
                    "-s",
                    str(FIXTURES / "tied_hubs.txt"),
                    "-S",
                    "17",
                    "-t",
                    "0",
                    "-e",
                    "1",
                    *extra_args,
                    "-o",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

            actual = {}
            for line in output.read_text(encoding="utf-8").splitlines():
                if not line or line.startswith(">"):
                    continue
                fields = line.split("\t")
                source = fields[0]
                self.assertEqual(len(fields[1:]) % 2, 0)
                for index in range(1, len(fields), 2):
                    actual[(source, fields[index])] = float(fields[index + 1])

            expected = {}
            with (FIXTURES / reference_name).open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    expected[(row["source"], row["target"])] = float(row["mi"])

            self.assertEqual(set(actual), set(expected))
            for edge, expected_mi in expected.items():
                self.assertAlmostEqual(actual[edge], expected_mi, places=6, msg=edge)

    def test_tied_bootstrap_seed_17_matches_reference_mi(self):
        self.assert_fixture_matches_reference(
            "tied_seed17_reference.tsv", "-r", "1"
        )

    def test_conditional_bootstrap_ranks_match_reference_mi(self):
        self.assert_fixture_matches_reference(
            "tied_conditional_bootstrap_seed17_reference.tsv",
            "-c",
            "-_FOXP1",
            "0.75",
            "-r",
            "1",
        )


if __name__ == "__main__":
    unittest.main()
