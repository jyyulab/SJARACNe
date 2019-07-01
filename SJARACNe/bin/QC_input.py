#!/usr/bin/env python3

import argparse
import sys

def main():
    head_description = 'Validating input files\n'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-e', '--exp-file', metavar='STR', required=True, help='exp file')
    parser.add_argument('-g', '--probe-file', metavar='STR', required=True, help='probe file')
    args = parser.parse_args()
    check_exp(args.exp_file)
    check_probe(args.probe_file)

def check_exp(input_file):
    with open(input_file, 'r') as fin:
        #process header
        header = fin.readline()
        words = header.split('\t')
        if words[0] != 'isoformId' or words[1] != 'geneSymbol':
            sys.exit('Error - Improper header in input file')
        #process rest of the file, making sure tabs are splitting entries
        for line in fin:
            words = line.split('\t')
            for word in words:
                if ' ' in word:
                    sys.exit('Error - spaces are not allowed, only tabs can delimit input file')
                if word.count('.') > 1:
                    sys.exit('Error - There are some numeric entries missing tab-spacing')

def check_probe(input_file):
    with open(input_file, 'r') as fin:
        #Make sure there is only one entry/gene per line
        for line in fin:
            size = len(line.split(' '))
            if size > 1:
                sys.exit('Error - There are more than one word per line in this probe file')


if __name__ == '__main__':
    main()
