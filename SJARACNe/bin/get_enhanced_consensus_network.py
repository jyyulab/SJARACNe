#!/usr/bin/env python3

import os
import sys
import re
from scipy import stats
import xlsxwriter
import pandas as pd
import pathlib
import argparse


def main():
    """ Handles arguments and invokes the driver function. """
    head_description = '''Add more information to a consensus network to create an enhanced network.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)
    parser.add_argument('-e', '--exp-mat', metavar='STR', required=True, help='expression matrix file')
    parser.add_argument('-n', '--network', metavar='STR', required=True, help='consensus network with 3 columns: '
                                                                              'isoformId1, isoformId2, '
                                                                              'mutual_information')
    parser.add_argument('-s', '--subnet', metavar='STR', help='gene symbols of interest')
    parser.add_argument('-o', '--out-dir', metavar='STR', required=True, help='output directory')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    create_enhanced_consensus_network(args.exp_mat, args.network, args.out_dir, args.subnet)
    print('All done', file=sys.stderr)


def create_enhanced_consensus_network(exp_mat, network, out_dir, subnet):
    """ Add more information to a consensus network to create an enhanced network.
    Args:
        exp_mat (str): path to an expression matrix file
        network (str): path to a consensus network file
        out_dir (str): path to an output directory
        subnet (str): path to a gene symbol file
    Return:
        None
    """
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # Build output file name based on input network file path
    network_file_name = os.path.basename(network)
    input_net_name_tokens = re.split("_|\.", network_file_name)
    out_file_name = ("_".join(input_net_name_tokens[0:2]) + "_ncol_" +
                     "_".join(input_net_name_tokens[3:len(input_net_name_tokens) - 1]))
    header = ("source", "target", "source.symbol", "target.symbol", "MI", "pearson",
              "spearman", "rho", "p-value")

    # Read in subnetwork file
    subnet_list = []
    out_subnet = None
    if subnet:
        out_subnet = open(pathlib.PurePath(out_dir).joinpath("sub_" + out_file_name + ".txt"), "w")
        with open(subnet) as subnet_file:
            for _id in subnet_file:
                subnet_list.append(_id.split("\n")[0].strip())

    exp = pd.read_csv(exp_mat, sep="\t", index_col=0)
    out_xlsx = xlsxwriter.Workbook(pathlib.PurePath(out_dir).joinpath(out_file_name + ".xlsx"))
    out_xlsx_sheet = out_xlsx.add_worksheet()
    row_index = 0
    out_xlsx_sheet.write_row(row_index, 0, header)

    with open(network, 'r') as fnet:
        fnet.readline()
        with open(pathlib.PurePath(out_dir).joinpath(out_file_name + ".txt"), 'w') as fout:
            fout.write('\t'.join(header) + '\n')
            for line in fnet:
                tokens = line.split('\t')
                node1 = tokens[0]
                node2 = tokens[1]
                mi = float(tokens[2])
                exp_values1 = exp.loc[node1].values[1:].astype(float)
                exp_values2 = exp.loc[node2].values[1:].astype(float)
                gene_symbol1 = exp.loc[node1].values[0]
                gene_symbol2 = exp.loc[node1].values[1]
                rho, intercept, r, p, stderr = stats.linregress(exp_values1, exp_values2)
                scc, sp = stats.spearmanr(exp_values1, exp_values2)
                pcc, pp = stats.pearsonr(exp_values1, exp_values2)

                row = (node1, node2, str(exp.loc[node1].values[0]), str(exp.loc[node2].values[0]),
                       "{0:.4f}".format(mi), "{0:.4f}".format(pcc), "{0:.4f}".format(scc),
                       "{0:.4f}".format(rho), "{0:.4f}".format(p))
                fout.write('\t'.join(row) + '\n')

                out_xlsx_sheet.write_row(row_index, 0, row)
                row_index += 1

                if out_subnet is not None:
                    if gene_symbol1 in subnet_list or gene_symbol2 in subnet_list:
                        out_subnet.write('\t'.join(row) + '\n')
    out_xlsx.close()

    if out_subnet is not None:
        out_subnet.close()
    print('All done', file=sys.stderr)


if __name__ == '__main__':
    main()
