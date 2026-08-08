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
            "C\tC\t1\t3\t2\t6\t4\t5\n"
            "D\tD\t2\t4\t1\t5\t3\t6\n"
            "E\tE\t3\t1\t5\t2\t6\t4\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.folder.cleanup()

    def run_sjaracne(self, adjacency_body, hubs=None, epsilon="0", tf_genes=None):
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
            epsilon,
            "-o",
            str(output),
        ]

        if hubs is not None:
            hub_file = self.workdir / "hubs.txt"
            hub_file.write_text(hubs, encoding="utf-8")
            command.extend(["-s", str(hub_file)])

        if tf_genes is not None:
            tf_file = self.workdir / "tf_genes.txt"
            tf_file.write_text(tf_genes, encoding="utf-8")
            command.extend(["-l", str(tf_file)])

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

    def test_sparse_rows_without_a_triangle_remain_unchanged(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.8\tD\t0.7\n"
            "C\n",
            hubs="A\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("[NETWORK] Applying DPI", result.stdout)
        self.assertEqual(
            self.data_rows(output),
            ["A\tB\t0.2\tC\t0.8\tD\t0.7"],
        )

    def test_selected_multi_hub_direct_intersection_prunes_triangle(self):
        # A-B has three stronger A neighbors but B has only two edges, so DPI
        # must scan B's direct row rather than the stronger-neighbor prefix.
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.9\tD\t0.8\tE\t0.7\n"
            "B\tC\t0.6\tD\t0.5\n",
            hubs="A\nB\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(output),
            [
                "A\tC\t0.9\tD\t0.8\tE\t0.7",
                "B\tC\t0.6\tD\t0.5",
            ],
        )

    def test_symmetric_all_gene_triangle_is_pruned_in_both_rows(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.8\n"
            "B\tA\t0.2\tC\t0.7\n"
            "C\tA\t0.8\tB\t0.7\n"
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(output),
            ["A\tC\t0.8", "B\tC\t0.7", "C\tA\t0.8\tB\t0.7"],
        )

    def test_sparse_symmetric_all_gene_intersection_prunes_weakest_edge(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.8\tD\t0.7\tE\t0.6\n"
            "B\tA\t0.2\tC\t0.5\n"
            "C\tA\t0.8\tB\t0.5\n"
            "D\tA\t0.7\n"
            "E\tA\t0.6\n"
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(output),
            [
                "A\tC\t0.8\tD\t0.7\tE\t0.6",
                "B\tC\t0.5",
                "C\tA\t0.8\tB\t0.5",
                "D\tA\t0.7",
                "E\tA\t0.6",
            ],
        )

    def test_tf_logic_protects_edge_from_non_tf_intermediate(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.8\n"
            "B\tC\t0.7\n",
            hubs="A\n",
            tf_genes="A\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.data_rows(output), ["A\tB\t0.2\tC\t0.8"])

    def test_tf_logic_also_protects_candidate_endpoint(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.9\tD\t0.8\n"
            "B\tC\t0.7\n",
            hubs="A\n",
            tf_genes="B\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(output),
            ["A\tB\t0.2\tC\t0.9\tD\t0.8"],
        )

    def test_protected_strongest_intermediate_does_not_hide_later_tf(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.9\tD\t0.8\tE\t0.7\n"
            "B\tC\t0.7\tD\t0.7\n",
            hubs="A\n",
            tf_genes="A\nD\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(output),
            ["A\tC\t0.9\tD\t0.8\tE\t0.7"],
        )

    def test_dpi_tolerance_uses_strict_threshold(self):
        retained, retained_output = self.run_sjaracne(
            "A\tB\t0.25\tC\t0.75\tD\t0.7\n"
            "B\tC\t0.5\n",
            hubs="A\n",
            epsilon="0.5",
        )
        self.assertEqual(retained.returncode, 0, retained.stdout + retained.stderr)
        self.assertEqual(
            self.data_rows(retained_output),
            ["A\tB\t0.25\tC\t0.75\tD\t0.7"],
        )

        pruned, pruned_output = self.run_sjaracne(
            "A\tB\t0.25\tC\t0.75\tD\t0.7\n"
            "B\tC\t0.500001\n",
            hubs="A\n",
            epsilon="0.5",
        )
        self.assertEqual(pruned.returncode, 0, pruned.stdout + pruned.stderr)
        self.assertEqual(
            self.data_rows(pruned_output),
            ["A\tC\t0.75\tD\t0.7"],
        )

        equal_ac, equal_ac_output = self.run_sjaracne(
            "A\tB\t0.25\tC\t0.5\tD\t0.8\tE\t0.7\n"
            "B\tC\t0.9\n",
            hubs="A\n",
            epsilon="0.5",
        )
        self.assertEqual(equal_ac.returncode, 0, equal_ac.stdout + equal_ac.stderr)
        self.assertEqual(
            self.data_rows(equal_ac_output),
            ["A\tB\t0.25\tC\t0.5\tD\t0.8\tE\t0.7"],
        )

    def test_nonempty_direct_row_prevents_reverse_fallback(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.9\tD\t0.8\tE\t0.7\n"
            "B\tA\t0.9\n"
            "C\tB\t0.6\n",
            hubs="A\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(output),
            ["A\tB\t0.2\tC\t0.9\tD\t0.8\tE\t0.7"],
        )

    def test_explicit_empty_row_uses_reverse_intersection_fallback(self):
        # B's effective incoming degree is two (A-B and C-B), smaller than the
        # three stronger A neighbors; this forces the reverse-index branch.
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.9\tD\t0.8\tE\t0.7\n"
            "B\n"
            "C\tB\t0.6\n",
            hubs="A\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            self.data_rows(output),
            ["A\tC\t0.9\tD\t0.8\tE\t0.7"],
        )

    def test_previously_pruned_edges_still_support_later_dpi(self):
        result, output = self.run_sjaracne(
            "A\tB\t0.2\tC\t0.4\tD\t0.9\tE\t0.8\n"
            "C\tB\t0.5\tD\t0.8\n",
            hubs="A\n",
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.data_rows(output), ["A\tD\t0.9\tE\t0.8"])

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
