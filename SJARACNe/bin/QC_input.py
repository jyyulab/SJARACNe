#!/usr/bin/env python3

import argparse
import sys
import logging


def main():
    head_description = 'Validating input files\n'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-e', '--exp-file', metavar='STR', required=True, help='exp file')
    parser.add_argument('-g', '--probe-file', metavar='STR', required=True, help='probe file')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    check_exp(args.exp_file)
    check_probe(args.probe_file)


def check_exp(input_file):
    total_genes = 0
    with open(input_file, 'r') as fin:
        # process header
        header = fin.readline()
        words = header.split('\t')
        if words[0] != 'isoformId' or words[1] != 'geneSymbol':
            sys.exit('Error - Improper header in input file: first two column names must be isoformId and '
                     'geneSymbol respectively.')
        entries_per_line = len(words)
        # process rest of the file, making sure tabs are splitting entries
        for line in fin:
            words = line.split('\t')
            total_genes += 1  # add one gene to the total count per line of the exp file
            if len(words) != entries_per_line:
                logging.info("Line {} does not have an appropriate number of entries".format(total_genes+1))
                sys.exit('Error - number of entries per line is not consistent across file. See line {}'.format(
                    total_genes+1))
            for word in words[2:]:
                if ' ' in word:
                    logging.info("Word with spaces is: {}".format(word))
                    sys.exit('Error - spaces are not allowed, only tabs can delimit input file. Space '
                             'found in line {}'.format(total_genes+1))
                if word.count('.') > 1 and word.isnumeric():
                    logging.info("Numeric entry missing spacing is: {}".format(word))
                    sys.exit('Error - There are some numeric entries missing tab-spacing in line '
                             '{}'.format(total_genes+1))
    logging.info("Number of genes in expression matrix: {}".format(total_genes))


def check_probe(input_file):
    probe_size = 0
    with open(input_file, 'r') as fin:
        # Make sure there is only one entry/gene per line
        for line in fin:
            probe_size += 1  # add one gene to size of probe file
            size = len(line.split(' '))
            if size > 1:
                sys.exit('Error - There are more than one word per line in this probe file')
    logging.info("Number of hub genes in probe file: {}".format(probe_size))


if __name__ == '__main__':
    main()

