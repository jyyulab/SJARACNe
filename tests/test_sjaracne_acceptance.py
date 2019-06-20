import unittest
import tempfile
import filecmp
import shlex
import subprocess

class TestSJARACNe(unittest.TestCase):
    def test_acceptance(self):
        cmd = 'sjaracne local -e ./tests/inputs/Tcell1170.exp -g ./tests/inputs/TcellTF.txt -o ./tests/test_run'
        exe = shlex.split(cmd)
        subprocess.check_call(exe)
        self.assertTrue(filecmp.cmp('./tests/answerkey/acceptance/consensus_network_ncol_.txt', './tests/test_run/consensus_network_ncol_.txt'))

if __name__ == '__main__':
    unittest.main()

