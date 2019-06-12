#!/usr/bin/env python3

import unittest
import filecmp
import os
from SJARACNe.bin.ch_line_ending import ch_line_ending as ch


class TestLineEnding(unittest.TestCase):
    def setUp(self):
        #create temporary input files that will hold line endings for their respective platforms
        self.file_names = ['win', 'mac', 'unix', 'inv', 'answer', 'out']
        self.endings = [b'\r\n', b'\r', b'\n', b'', b'\n']
        for i in range(5):
            with open(self.file_names[i], 'wb') as fin:
                fin.write(b'Hello world!' + self.endings[i])
                fin.write(b'This is an example file')
                fin.close()
        
        #create empty output file
        with open('out', 'wb') as fin:
            fin.write(b'')
            fin.close()

    def test_infile_same_as_outfile(self):
        with self.assertRaises(SystemExit) as err:
            ch('win', 'win')
        self.assertEqual(err.exception.code, 'Error - you must omit output file argument if it is identical to input file')
        

    def test_first_line_is_invalid(self):
        with self.assertRaises(SystemExit) as err:
            ch('inv', 'out')

        with open('inv', 'rb') as fin:
            first_line = fin.readline()
            fin.close();
        self.assertEqual(err.exception.code, 'Error - invalid line ending in the first line: {}'.format(first_line))


    def test_unix_ending_file(self):
        output_file = ch('unix', 'out')      
        self.assertTrue(filecmp.cmp(output_file, 'answer'))

    def test_windows_ending_file(self):
        output_file = ch('win', 'out')
        self.assertTrue(filecmp.cmp(output_file, 'answer'))


    def test_mac_ending(self):
        output_file = ch('mac', 'out')
        self.assertTrue(filecmp.cmp(output_file, 'answer'))

    def tearDown(self):
        #deletes all files that have been used
        for name in self.file_names:
            os.remove(name)


if __name__ == "__main__":
    unittest.main()
