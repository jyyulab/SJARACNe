#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import shlex
import logging
import pathlib


def main():
    head_description = '''SJARACNe is a scalable tool for gene network reverse engineering.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)

    # Create a parent parser with common arguments for every subparser
    parent_parser = argparse.ArgumentParser(description='hello', add_help=False)
    parent_parser.add_argument('-e', '--exp-file', metavar='FILE', required=True,
                               help='Path to an expression matrix file, row indexes are used as the'
                               'nodes in the network.')
    parent_parser.add_argument('-g', '--hub-genes', metavar='FILE', required=True,
                               help='Path to a file containing a list of symbols to be considered as hub genes.')
    parent_parser.add_argument('-pc', '--p-value-consensus', metavar='FLOAT', default=1e-5,
                               help='P-value threshold to select edges in building consensus network.')
    parent_parser.add_argument('-pb', '--p-value-bootstrap', metavar='FLOAT', default=1e-7,
                               help='P-value threshold to filter mutual information in building bootstrap networks.')
    parent_parser.add_argument('-d', '--depth', metavar='INT', default=40, help='maximum partitioning depth.')
    parent_parser.add_argument('-c', '--config-dir', metavar='DIR', help='Directory containing ARACNe configuration '
                                                                         'files. Use default configs if not provided.')
    parent_parser.add_argument('-n', '--bootstrap-num', metavar='INT', default=100,
                               help='Number of bootstrap networks to generate.')
    parent_parser.add_argument('-o', '--output-dir', metavar='DIR', required=True,
                               help='Path to final output directory.')

    subparsers = parser.add_subparsers(title='Subcommands', help='platforms', dest='subcommand')
    subparsers.required = True
    # Create a subparser for running cwltool
    subparser_local = subparsers.add_parser('local', parents=[parent_parser], help='run cwltool in a local workstation')
    subparser_local.add_argument('-s', '--serial', help='run cwltool in serial mode', action='store_true')

    # Create a subparser for running cwlexec
    subparser_lsf = subparsers.add_parser('lsf', parents=[parent_parser], help='run cwlexec in a IBM LSF cluster')
    subparser_lsf.add_argument('-j', '--config-json', metavar='FILE', required=True, help='LSF-specific configuration '
                                                                                         'file in JSON format to be '
                                                                                         'used for workflow execution')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    # to make executable and config findable
    installed_path = os.path.dirname(os.path.realpath(__file__))
    os.environ['PATH'] += (os.pathsep + installed_path + '/bin')
    cwl_path = installed_path + '/cwl'
    default_config_path = installed_path + '/config'
    if args.config_dir:
        config_dir = args.config_dir
    else:
        config_dir = default_config_path

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    output_dir_name = os.path.basename(args.output_dir)
    # Create input yml file in a temp directory
    with open(pathlib.PurePath(args.output_dir).joinpath('sjaracne_workflow.yml'), 'w') as fp_yml:
        logging.info(fp_yml.name)
        contents = 'exp_file:\n  class: File\n  path: {}\n' \
                   'probe_file:\n  class: File\n  path: {}\n' \
                   'p_value_consensus: {}\n' \
                   'p_value_bootstrap: {}\n' \
                   'depth: {}\n' \
                   'aracne_config_dir:\n  class: Directory\n  path: {}\n' \
                   'bootstrap_num: {}\n' \
                   'final_out_dir_name: {}'.format(os.path.abspath(args.exp_file), os.path.abspath(args.hub_genes),
                                                   args.p_value_consensus, args.p_value_bootstrap, args.depth,
                                                   config_dir, args.bootstrap_num, output_dir_name)
        logging.info(contents)
        fp_yml.write(contents)
        fp_yml.flush()
        fp_yml.seek(0)

        if args.subcommand == 'local':
            if args.serial:
                cmd = 'cwltool --outdir {} {}/sjaracne_workflow.cwl {}'.format(args.output_dir, cwl_path,
                                                                               fp_yml.name)
            else:
                cmd = 'cwltool --parallel --outdir {} {}/sjaracne_workflow.cwl {}'.format(args.output_dir,
                                                                                          cwl_path, fp_yml.name)
        elif args.subcommand == 'lsf':
                cmd = 'cwlexec -pe PATH -c {} --outdir {} {}/sjaracne_workflow.cwl {}'.format(
                    args.config_json, args.output_dir, cwl_path, fp_yml.name)
        else:
            sys.exit('Error - invalid subcommand.')
        logging.info(cmd)
        run_shell_command_call(cmd)

    logging.info('All done.')


def run_shell_command_call(cmd):
    """ Wrapper of subprocess.check_call to take a cmd string as input
    Args:
        cmd (str): command to run
    """
    cmd_to_exec = shlex.split(cmd)
    subprocess.check_call(cmd_to_exec)


if __name__ == "__main__":
    main()
