#!/usr/bin/env python3

import unittest
import filecmp
import tempfile
import sys
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

if __name__ == '__main__':
    unittest.main()
