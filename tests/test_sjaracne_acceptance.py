#!/usr/bin/env python3

import unittest
import tempfile
import filecmp
import shlex
import subprocess


class TestSJARACNe(unittest.TestCase):
    def test_acceptance(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = 'sjaracne local -e ./tests/inputs/Tcell1170.exp -g ./tests/inputs/TcellTF.txt -n 5 -o {}'.format(tmpdir)
            exe = shlex.split(cmd)
            subprocess.check_call(exe)
            cmd2 = 'diff ./tests/answerkey/acceptance/cnn_5.txt {}/consensus_network_ncol_.txt'.format(tmpdir)
            exe2 = shlex.split(cmd2)
            subprocess.check_call(exe2)
            self.assertTrue(filecmp.cmp('./tests/answerkey/acceptance/cnn_5.txt', '{}/consensus_network_ncol_.txt'.format(tmpdir)))


if __name__ == '__main__':
    unittest.main()

