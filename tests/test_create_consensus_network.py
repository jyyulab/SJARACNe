#!/usr/bin/env python3

import unittest
import filecmp
import tempfile
import sys
from pathlib import Path
from SJARACNe.bin.create_consensus_network import create_consensus_network as cn
from SJARACNe.bin.create_consensus_network import create_enhanced_consensus_network as ecn
from SJARACNe.bin.create_consensus_network import uprob


class TestConsensusNetwork(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.folder = tempfile.TemporaryDirectory()
        cn('./tests/inputs/adjmat_dir', 0.05, self.folder.name)

    def test_consensus_network_3col(self):
        self.assertTrue(filecmp.cmp('./tests/answerkey/consensus_network_3col_.txt', self.folder.name + '/consensus_network_3col_.txt'))
        
    def test_bootstrap_info(self):
        self.assertTrue(filecmp.cmp('./tests/answerkey/bootstrap_info_.txt', self.folder.name + '/bootstrap_info_.txt'))

    def test_parameter_info(self):  
        #delete last line of parameter info, as it will not match with answer key
        with open(self.folder.name + '/parameter_info_.txt', 'r') as info:
            lines = info.readlines()
            lines = lines[:-1]
            info.close()
        with open(self.folder.name + '/parameter_info_.txt', 'w') as info:
            info.writelines(lines)
            info.close()
        self.assertTrue(filecmp.cmp('./tests/answerkey/parameter_info_.txt', self.folder.name + '/parameter_info_.txt'))

    '''    
    def test_enhanced_consensus_network(self):
        ecn('./inputs/BRCA100.exp', './answerkey/consensus_network_3col_.txt', self.folder.name + '/consensus_network_ncol_.txt')
        self.assertTrue(filecmp.cmp('./answerkey/consensus_network_ncol_.txt', self.folder.name + '/consensus_network_ncol_.txt'))
    '''

    #Testing different values on the function uprob(z-score)
    def test_uprob_100(self):
        self.assertAlmostEqual(0, uprob(100))

    def test_uprob_2(self):
        self.assertAlmostEqual(0.022750150253, uprob(2))

    def test_uprob_half(self):
        self.assertAlmostEqual(0.308518356555, uprob(0.5))

    def test_uprob_neg1(self):
        self.assertAlmostEqual(0.841344680778, uprob(-1))


class TestConsensusDuplicateEdges(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.workdir = Path(self.folder.name)

    def tearDown(self):
        self.folder.cleanup()

    @staticmethod
    def write_bootstrap(directory, filename, body):
        directory.mkdir(parents=True, exist_ok=True)
        (directory / filename).write_text(
            ">  Input file test.exp\n" + body,
            encoding="utf-8",
        )

    @staticmethod
    def read_network(output_dir):
        return (output_dir / "consensus_network_3col_.txt").read_text(
            encoding="utf-8"
        )

    def test_duplicate_edge_occurrences_count_once_per_bootstrap(self):
        clean_dir = self.workdir / "clean"
        duplicate_dir = self.workdir / "duplicate"
        clean_output = self.workdir / "clean_output"
        duplicate_output = self.workdir / "duplicate_output"

        self.write_bootstrap(clean_dir, "run_1.adj", "A\tB\t0.5\tC\t0.2\n")
        self.write_bootstrap(clean_dir, "run_2.adj", "A\tB\t0.7\n")
        self.write_bootstrap(
            duplicate_dir,
            "run_1.adj",
            "A\tB\t0.5\tB\t0.5\tC\t0.2\nA\tB\t0.5\n",
        )
        self.write_bootstrap(duplicate_dir, "run_2.adj", "A\tB\t0.7\n")

        cn(clean_dir, 1.0, clean_output)
        with self.assertLogs(level="WARNING") as logs:
            cn(duplicate_dir, 1.0, duplicate_output)

        self.assertIn("Ignored 2 duplicate edge occurrence(s)", "\n".join(logs.output))
        self.assertEqual(
            self.read_network(duplicate_output), self.read_network(clean_output)
        )
        self.assertIn("A\tB\t0.6000\n", self.read_network(duplicate_output))
        self.assertEqual(
            (duplicate_output / "bootstrap_info_.txt").read_text(encoding="utf-8"),
            (clean_output / "bootstrap_info_.txt").read_text(encoding="utf-8"),
        )

    def test_conflicting_duplicate_mi_values_fail_closed(self):
        adjacency_dir = self.workdir / "conflicting"
        output_dir = self.workdir / "output"
        self.write_bootstrap(
            adjacency_dir,
            "run_1.adj",
            "A\tB\t0.5\nA\tB\t0.6\n",
        )

        with self.assertRaisesRegex(
            ValueError,
            "Conflicting MI values for duplicate edge A----B.*0.5 versus 0.6",
        ):
            cn(adjacency_dir, 1.0, output_dir)

        self.assertFalse((output_dir / "consensus_network_3col_.txt").exists())

    def test_equivalent_mi_text_is_an_exact_numeric_duplicate(self):
        adjacency_dir = self.workdir / "equivalent"
        output_dir = self.workdir / "output"
        self.write_bootstrap(
            adjacency_dir,
            "run_1.adj",
            "A\tB\t0.5\nA\tB\t0.5000\n",
        )

        with self.assertLogs(level="WARNING"):
            cn(adjacency_dir, 1.0, output_dir)

        self.assertIn("A\tB\t0.5000\n", self.read_network(output_dir))

    def test_nonfinite_mi_fails_closed(self):
        adjacency_dir = self.workdir / "nonfinite"
        output_dir = self.workdir / "output"
        self.write_bootstrap(adjacency_dir, "run_1.adj", "A\tB\tnan\n")

        with self.assertRaisesRegex(
            ValueError,
            "Non-finite MI value for edge A----B.*nan",
        ):
            cn(adjacency_dir, 1.0, output_dir)

        self.assertFalse((output_dir / "consensus_network_3col_.txt").exists())

    def test_reverse_ordered_edge_is_not_treated_as_a_duplicate(self):
        adjacency_dir = self.workdir / "directed"
        output_dir = self.workdir / "output"
        self.write_bootstrap(
            adjacency_dir,
            "run_1.adj",
            "A\tB\t0.5\nB\tA\t0.7\n",
        )

        cn(adjacency_dir, 1.0, output_dir)

        network = self.read_network(output_dir)
        self.assertIn("A\tB\t0.5000\n", network)
        self.assertIn("B\tA\t0.7000\n", network)

if __name__ == '__main__':
    unittest.main()
