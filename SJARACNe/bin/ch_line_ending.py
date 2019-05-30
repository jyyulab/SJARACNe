#!/usr/bin/env python3

import os
import sys
import argparse


def main():
    """ Handles arguments and invokes the driver function. """
    head_description = 'Convert line ending to \n'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-i', '--input-file', metavar='STR', required=True, help='Input file')
    parser.add_argument('-o', '--output-file', metavar='STR', help='Output file')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    ch_line_ending(args.input_file, args.output_file)


def ch_line_ending(input_file, output_file):
    """ Replace windows line ending '\r\n' or Mac line ending '\r' with Unix line ending '\n'
        Args:
            input_file (str): input file path
            output_file (str or None): output file path. Same as input_file if it is None
        Returns:
            output_file (str): file with Unix line ending
    """
    if output_file == input_file:
        sys.exit('Error - you must omit output file argument if it is identical to input file')

    windows_line_ending = b'\r\n'
    mac_line_ending = b'\r'
    unix_line_ending = b'\n'

    with open(input_file, 'rb') as fin:
        first_line = fin.readline()
        if windows_line_ending in first_line:
            print('Windows line ending', file=sys.stderr)
        elif mac_line_ending in first_line:
            print('Mac line ending', file=sys.stderr)
        elif unix_line_ending in first_line:
            print('Unix line ending. No need to convert', file=sys.stderr)
            fin.close()
            return input_file
        else:
            sys.exit('Error - invalid line ending in the first line: {}'.format(first_line))

    if output_file is None:
        output_file = '{}.tmp'.format(input_file)
    with open(output_file, 'w', newline='\n') as fout:
        with open(input_file, 'r') as fin:
            for line in fin:
                fout.write(line)

    if output_file is None:
        os.remove(input_file)
        os.rename(output_file, input_file)
        output_file = input_file
    print('Done', file=sys.stderr)
    return output_file


if __name__ == '__main__':
    main()
