#!/usr/bin/env python3

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


def find_sjaracne_executable():
    configured = os.environ.get("SJARACNE_TEST_EXE")
    candidates = [
        configured,
        Path(__file__).resolve().parents[1] / "SJARACNe" / "bin" / "sjaracne.exe",
        shutil.which("sjaracne.exe"),
    ]

    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(candidate)
    return None


SJARACNE_EXE = find_sjaracne_executable()


@unittest.skipUnless(SJARACNE_EXE, "sjaracne.exe is not built")
class TestAdjacencyInput(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)
        self.expression = self.workdir / "small.exp"
        self.expression.write_text(
            "isoformId\tgeneSymbol\ts1\ts2\ts3\ts4\ts5\ts6\n"
            "A\tA\t1\t2\t3\t4\t5\t6\n"
            "B\tB\t6\t5\t4\t3\t2\t1\n"
            "C\tC\t1\t3\t2\t6\t4\t5\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.folder.cleanup()

    def run_sjaracne(self, adjacency_body, hubs=None):
        adjacency = self.workdir / "input.adj"
        adjacency.write_text(adjacency_body, encoding="utf-8")
        output = self.workdir / "output.adj"
        command = [
            SJARACNE_EXE,
            "-i",
            str(self.expression),
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
            hub_file = self.workdir / "hubs.txt"
            hub_file.write_text(hubs, encoding="utf-8")
            command.extend(["-s", str(hub_file)])

        result = subprocess.run(command, capture_output=True, text=True, check=False)
        return result, output

    @staticmethod
    def data_rows(output):
        return [
            line
            for line in output.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith(">")
        ]

    def test_requested_hub_missing_as_source_row_fails_cleanly(self):
        result, output = self.run_sjaracne("A\tB\t0.5\n", hubs="C\n")

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(
            "does not contain a source row for requested hub: C", result.stderr
        )
        self.assertIn(
            "Refusing to treat absent adjacency rows as empty networks",
            result.stderr,
        )
        self.assertFalse(output.exists())

    def test_requested_hub_appearing_only_as_target_is_still_absent(self):
        result, output = self.run_sjaracne("A\tC\t0.5\n", hubs="C\n")

        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(
            "does not contain a source row for requested hub: C", result.stderr
        )
        self.assertFalse(output.exists())

    def test_explicit_source_only_row_is_a_valid_empty_row(self):
        result, output = self.run_sjaracne("A\tB\t0.5\nC\n", hubs="C\n")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(output.is_file())
        self.assertEqual(self.data_rows(output), [])

    def test_sparse_adjacency_is_safe_in_all_gene_mode(self):
        result, output = self.run_sjaracne("A\tB\t0.5\n")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(output.is_file())
        self.assertEqual(self.data_rows(output), ["A\tB\t0.5"])

    def test_three_gene_dpi_triangle_prunes_the_weakest_edge(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.8\n"
            "B\tC\t0.7\n"
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[NETWORK] Applying DPI", result.stdout)
        self.assertTrue(output.is_file())
        self.assertEqual(
            self.data_rows(output),
            ["A\tC\t0.8", "B\tC\t0.7"],
        )

    def test_single_imported_source_row_skips_impossible_dpi(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.8\n",
            hubs="A\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[NETWORK] Skipping DPI", result.stdout)
        self.assertNotIn("[NETWORK] Applying DPI", result.stdout)
        self.assertEqual(self.data_rows(output), ["A\tB\t0.2\tC\t0.8"])

    def test_imported_extra_source_row_keeps_dpi_available(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.8\n"
            "B\tC\t0.7\n",
            hubs="A\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[NETWORK] Applying DPI", result.stdout)
        self.assertNotIn("[NETWORK] Skipping DPI", result.stdout)
        self.assertEqual(self.data_rows(output), ["A\tC\t0.8"])

    def test_duplicate_accession_resolves_to_first_expression_row(self):
        self.expression.write_text(
            "isoformId\tgeneSymbol\ts1\ts2\ts3\ts4\ts5\ts6\n"
            "DUP\tDUP_FIRST\t1\t2\t3\t4\t5\t6\n"
            "Y\tY\t6\t5\t4\t3\t2\t1\n"
            "DUP\tDUP_SECOND\t1\t3\t2\t6\t4\t5\n"
            "X\tX\t2\t1\t4\t3\t6\t5\n",
            encoding="utf-8",
        )

        result, output = self.run_sjaracne("X\tDUP\t0.5\tY\t0.6\n")

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.data_rows(output), ["X\tDUP\t0.5\tY\t0.6"])


if __name__ == "__main__":
    unittest.main()
