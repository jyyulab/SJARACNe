#!/usr/bin/env python3

import unittest
import tempfile
import subprocess
from pathlib import Path


class TestSJARACNe(unittest.TestCase):
    def test_acceptance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir)
            output_dir = workdir / 'output'
            temp_dir = workdir / 'tmp'
            subprocess.check_call([
                'sjaracne',
                'local',
                '-e',
                './tests/inputs/Tcell1170.exp',
                '-g',
                './tests/inputs/TcellTF.txt',
                '-n',
                '5',
                '-pc',
                '0.01',
                '-o',
                str(output_dir),
                '-tmp',
                str(temp_dir),
            ])

            expected = Path('./tests/answerkey/acceptance/cnn_5.txt')
            actual = output_dir / 'consensus_network_ncol_.txt'
            self.assertEqual(
                expected.read_text(encoding='utf-8'),
                actual.read_text(encoding='utf-8'),
            )


if __name__ == '__main__':
    unittest.main()

