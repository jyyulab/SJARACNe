#!/usr/bin/env python3
import unittest
import filecmp
import tempfile
import os
from SJARACNe.bin.ch_line_ending import ch_line_ending as ch


class TestLineEnding(unittest.TestCase):
    @classmethod
    def setUpClass(self): 
        #create empty output file
        self.out = tempfile.NamedTemporaryFile()
        self.out.write(b'')
        #create answer file
        self.answer = tempfile.NamedTemporaryFile()
        self.answer.write(b'Hello world!' + b'\n')
        self.answer.write(b'This is an example file')
        self.answer.seek(0)

    def test_infile_same_as_outfile(self):
        fp = tempfile.NamedTemporaryFile()
        fp.write(b'Hello world!\r\nThis is an example file\r\n')
        with self.assertRaises(SystemExit) as err:
            ch(fp.name, fp.name)
        self.assertEqual(err.exception.code, 'Error - you must omit output file argument if it is identical to input file')
        

    def test_first_line_is_invalid(self):
        fp = tempfile.NamedTemporaryFile(delete=False)
        fp.write(b'Hello world!')
        fp.write(b'This is an example file')
        with self.assertRaises(SystemExit) as err:
            ch(fp.name, self.out.name)

        with open(fp.name, 'rb') as fin:
            first_line = fin.readline()
        self.assertEqual(err.exception.code, 'Error - invalid line ending in the first line: {}'.format(first_line))


    def test_unix_ending_file(self):
        fp = tempfile.NamedTemporaryFile(delete=False)
        fp.write(b'Hello world!' + b'\n')
        fp.write(b'This is an example file')
        fp.seek(0)
        output_file = ch(fp.name, self.out.name)
        self.assertTrue(filecmp.cmp(output_file, self.answer.name))

    def test_windows_ending_file(self):
        fp = tempfile.NamedTemporaryFile(delete=False) 
        fp.write(b'Hello world!' + b'\r\n')
        fp.write(b'This is an example file')
        fp.seek(0)
        output_file = ch(fp.name, self.out.name)
        self.assertTrue(filecmp.cmp(output_file, self.answer.name))

    def test_mac_ending(self):
        fp = tempfile.NamedTemporaryFile(delete=False)
        fp.write(b'Hello world!' + b'\r')
        fp.write(b'This is an example file')
        fp.seek(0)
        output_file = ch(fp.name, self.out.name)
        self.assertTrue(filecmp.cmp(output_file, self.answer.name))
    
    @classmethod
    def tearDownClass(self):
        self.out.close()
        self.answer.close()

if __name__ == "__main__":
    unittest.main()
