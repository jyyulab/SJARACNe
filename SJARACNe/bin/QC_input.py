#!/usr/bin/env python3

import argparse
import logging
import re
import sys


DEFAULT_REPORT_FILE = 'hub_overlap_validation.txt'
BOUNDARY_WHITESPACE = ' \t\n\r\f\v'


def main():
    head_description = 'Validating input files\n'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=head_description)
    parser.add_argument('-e', '--exp-file', metavar='STR', required=True, help='exp file')
    parser.add_argument('-g', '--probe-file', metavar='STR', required=True, help='probe file')
    parser.add_argument('-o', '--output-file', metavar='STR', default=DEFAULT_REPORT_FILE,
                        help='output validation report')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    validate_inputs(args.exp_file, args.probe_file, args.output_file)


def check_exp(input_file):
    total_genes = 0
    expression_ids = set()

    with open(input_file, 'r', encoding='utf-8-sig', newline=None) as fin:
        # process header
        header = fin.readline().rstrip('\r\n')
        words = header.split('\t')
        if len(words) < 2 or words[0] != 'isoformId' or words[1] != 'geneSymbol':
            sys.exit('Error - Improper header in input file: first two column names must be isoformId and '
                     'geneSymbol respectively.')
        entries_per_line = len(words)
        # process rest of the file, making sure tabs are splitting entries
        for line in fin:
            line = line.rstrip('\r\n')
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
            expression_ids.add(words[0])
    logging.info("Number of genes in expression matrix: {}".format(total_genes))
    return expression_ids


def check_probe(input_file):
    with open(input_file, 'rb') as fin:
        try:
            contents = fin.read().decode('utf-8-sig')
        except UnicodeDecodeError as error:
            sys.exit('Error - probe file is not valid UTF-8: {}'.format(error))

    probe_ids = []
    seen = set()
    duplicate_count = 0

    # Match the C++ hub-list parser: accept LF, CRLF, lone CR, and mixed endings;
    # trim ASCII boundary whitespace; skip blank records; preserve case and
    # internal whitespace; and keep the first occurrence of each identifier.
    for line in re.split(r'\r\n|\r|\n', contents):
        probe_id = line.strip(BOUNDARY_WHITESPACE)
        if not probe_id:
            continue
        if probe_id in seen:
            duplicate_count += 1
            continue
        seen.add(probe_id)
        probe_ids.append(probe_id)

    logging.info("Number of unique hub genes in probe file: {}".format(len(probe_ids)))
    if duplicate_count:
        logging.info("Duplicate hub genes ignored: {}".format(duplicate_count))
    return probe_ids, duplicate_count


def validate_hub_overlap(expression_ids, probe_ids):
    matched = [probe_id for probe_id in probe_ids if probe_id in expression_ids]
    missing = [probe_id for probe_id in probe_ids if probe_id not in expression_ids]

    logging.info(
        "Hub overlap: requested={}, matched={}, missing={}".format(
            len(probe_ids), len(matched), len(missing)
        )
    )

    if missing:
        preview = ', '.join(missing[:20])
        suffix = '' if len(missing) <= 20 else ', ...'
        logging.warning("Hub genes absent from expression matrix: {}{}".format(preview, suffix))

    if not matched:
        sys.exit(
            'Error - zero hub genes from the probe file matched the expression matrix '
            '(requested: {}). Refusing to start bootstrap jobs.'.format(len(probe_ids))
        )

    return matched, missing


def write_validation_report(output_file, expression_ids, probe_ids, duplicate_count,
                            matched, missing):
    with open(output_file, 'w', encoding='utf-8', newline='\n') as fout:
        fout.write('status\tpassed\n')
        fout.write('expression_genes\t{}\n'.format(len(expression_ids)))
        fout.write('hub_genes_requested\t{}\n'.format(len(probe_ids)))
        fout.write('hub_genes_matched\t{}\n'.format(len(matched)))
        fout.write('hub_genes_missing\t{}\n'.format(len(missing)))
        fout.write('hub_duplicates_ignored\t{}\n'.format(duplicate_count))
        fout.write('matched_hubs\t{}\n'.format(','.join(matched)))
        fout.write('missing_hubs\t{}\n'.format(','.join(missing)))


def validate_inputs(exp_file, probe_file, output_file=DEFAULT_REPORT_FILE):
    expression_ids = check_exp(exp_file)
    probe_ids, duplicate_count = check_probe(probe_file)
    matched, missing = validate_hub_overlap(expression_ids, probe_ids)
    write_validation_report(
        output_file, expression_ids, probe_ids, duplicate_count, matched, missing
    )
    logging.info("Hub-overlap validation passed; report: {}".format(output_file))
    return matched, missing


if __name__ == '__main__':
    main()

