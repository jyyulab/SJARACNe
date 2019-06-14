#!/usr/env/bin python3

import unittest
import filecmp
import sys
import os
import shutil
from SJARACNe.bin.create_consensus_network import create_consensus_network as cn
from SJARACNe.bin.create_consensus_network import create_enhanced_consensus_network as ecn
from SJARACNe.bin.create_consensus_network import uprob

class TestConsensusNetwork(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.path = './test_answers'
        if self.path is None:
            os.mkdir('./test_answers')
        cn('./tests/inputs/adjmat_dir', 0.05, './test_answers')

    
    def test_consensus_network_3col(self):
        self.assertTrue(filecmp.cmp('./tests/answerkey/consensus_network_3col_.txt', './test_answers/consensus_network_3col_.txt'))
        
    def test_bootstrap_info(self):
        self.assertTrue(filecmp.cmp('./tests/answerkey/bootstrap_info_.txt', './test_answers/bootstrap_info_.txt'))

    def test_parameter_info(self):  
        #delete last line of parameter info, as it will not match with answer key
        with open('./tests/test_answers/parameter_info_.txt', 'r') as info:
            lines = info.readlines()
            lines = lines[:-1]
            info.close()
        with open('./test_answers/parameter_info_.txt', 'w') as info:
            info.writelines(lines)
            info.close()
        self.assertTrue(filecmp.cmp('./tests/answerkey/parameter_info_.txt', './test_answers/parameter_info_.txt'))


    '''
    def test_enhanced_consensus_network(self):
        ecn('./inputs/BRCA100.exp', './answerkey/consensus_network_3col_.txt', './test_answers/consensus_network_ncol_.txt')
        self.assertTrue(filecmp.cmp('./answerkey/consensus_network_ncol_.txt', './test_answers/consensus_network_ncol_.txt'))
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
    
    @classmethod
    def tearDownClass(self): 
        shutil.rmtree('./test_answers')
        
    

if __name__ == '__main__':
    unittest.main()
