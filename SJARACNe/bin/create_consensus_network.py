#!/usr/bin/env python3

import sys
import os
import argparse
import math
import logging
import numpy as np
import pathlib
import re
from scipy import stats
import pandas as pd


def main():
    """ Handles arguments and invokes the driver function. """
    head_description = '''Create a consensus network based on SJARACNe bootstrap networks.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=head_description)
    parser.add_argument('-a', '--adjmat-dir', metavar='STR', required=True, help='directory with adjacent matrix')
    parser.add_argument('-p', '--p-value', metavar='STR', required=True, help='P value threshold')
    parser.add_argument('-e', '--exp-mat', metavar='STR', required=True, help='expression matrix file')
    parser.add_argument('-o', '--out-dir', metavar='STR', required=True, help='output directory')
    parser.add_argument('-s', '--subnet', metavar='STR', help='file with gene symbols of interest to build a subnet')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    logging.basicConfig(level=logging.INFO)
    logging.info('Create an initial consensus network ...')
    network = create_consensus_network(args.adjmat_dir, args.p_value, args.out_dir)
    logging.info('Done')
    logging.info('Create an enhanced consensus network ...')
    create_enhanced_consensus_network(args.exp_mat, network, args.out_dir, args.subnet)
    logging.info('All done')


def create_consensus_network(adjmat_dir, p_value, out_dir):
    """ Create a consensus network based on SJARACNe bootstrap networks
    Args:
        adjmat_dir: directory with adjacent matrix
        p_value: P value threshold
        out_dir: output directory
    Returns:
        none
    """
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    total_edge_in_runs = []
    bootstrap_run_num = 0
    parameters = ''
    total_edge_number = {}
    total_mi = {}

    # Processing all bootstrap networks, summarizing them into corresponding variables
    for adj_file in os.listdir(adjmat_dir):
        total_edge_in_runs.append(0)
        # Opening each bootstrap file
        with open(pathlib.PurePath(adjmat_dir).joinpath(adj_file), "r") as fadj:
            for line in fadj:
                # Processing header lines
                if line[0] == '>' and bootstrap_run_num == 0:
                    parameters += line
                # Processing non header lines representing the network
                if line[0] != '>':
                    tokens = line.split('\t')
                    # Tokenizing each non header line with tab delimiter
                    hub_id = tokens[0]  # First token is the hub id
                    # Iterating on all adjacent genes: Odd indexes are the connected genes and even indexes are
                    # the corresponding value to the edge between the hub gene and the gene with an odd index
                    # appearing before the value in the tokens list
                    for index in range(1, len(tokens), 2):
                        key = hub_id + "----" + tokens[index]  # Creating a key for the edge
                        if key in total_edge_number:
                            # Updating the total number of edges observed for the particular key (edge)
                            total_edge_number[key] += 1
                            # Updating total MI between the genes involving in the particular key (edge)
                            total_mi[key] += float(tokens[index + 1])
                        else:
                            # Initializing the total number of edges observed for the particular key (edge)
                            # if the key (edge) is newly generated
                            total_edge_number[key] = 1
                            # Initializing total MI between the genes involving in the particular key (edge)
                            # if the key (edge) is newly generated
                            total_mi[key] = float(tokens[index + 1])
                        # Increment the total number of edges processed so far for a bootstrap run
                        total_edge_in_runs[bootstrap_run_num] += 1
        # Increment the bootstrap file index
        bootstrap_run_num += 1

    mu = 0
    sigma = 0
    # Computing mu and sigma across all bootstrap files
    for i in range(0, bootstrap_run_num):
        prob = float(total_edge_in_runs[i]) / float(len(total_edge_number))
        mu += prob
        sigma += prob * (1 - prob)
    sigma = np.sqrt(sigma)

    # Writing out the summary of all bootstrap files into bootstrap_info.txt file
    with open(pathlib.PurePath(out_dir).joinpath('bootstrap_info_.txt'), 'w') as f_info:
        f_info.write('Total edge tested: {}\n'.format(str(len(total_edge_number))))
        f_info.write('Bonferroni corrected (0.05) alpha: {}\n'.format(str(0.05 / len(total_edge_number))))
        f_info.write('mu: {}\n'.format(str(mu)))
        f_info.write('sigma: {}\n'.format(str(sigma)))

    # Setting p_threshold to the given value, if no given value, set to Bonferroni corrected value
    p_threshold = 0.05 / len(total_edge_number)
    if p_value is not None:
        p_threshold = float(p_value)

    # Writing out the parameters that the bootstrap networks are constructed with plus other
    # parameters that is used to create consensus network
    parameters += '>  Bootstrap No: {}\n'.format(str(bootstrap_run_num))
    parameters += '>  Source: sjaracne2\n'
    out_network_path = pathlib.PurePath(out_dir).joinpath('consensus_network_3col_.txt')
    parameters += '>  Output network: {}\n'.format(out_network_path)
    with open(pathlib.PurePath(out_dir).joinpath('parameter_info_.txt'), 'w') as parameter_file:
        parameter_file.write(parameters)

    # Writing out the consensus network preserving edges with statistically significant support
    with open(out_network_path, 'w') as f_consensus_network:
        header = 'source\ttarget\tMI\n'
        f_consensus_network.write(header)

        # Iterate over all edges in a sorted fashion
        for key in sorted(total_edge_number.keys()):
            # Extract first two gene involving an edge from the key (edge)
            tks = key.split('----')
            gene1 = tks[0]
            gene2 = tks[1]

            # Compute the z score of normal distribution
            z = float(total_edge_number[key] - mu) / float(sigma) if sigma != 0 else 100

            # Compute p-value corresponding to the z score
            pval = uprob(z)

            # Decision making if the edge has enough support or not and therefore if it has to be remained or removed
            if pval < p_threshold:
                # Computing MI corresponding to an edge remaining in the network
                mi = '{0:.4f}'.format(float(total_mi[key]) / float(total_edge_number[key]))
                f_consensus_network.write('{}\t{}\t{}\n'.format(gene1, gene2, mi))
    return out_network_path


def create_enhanced_consensus_network(exp_mat, network, out_dir, subnet=None):
    """ Add more information to a consensus network to create an enhanced network.
    Args:
        exp_mat (str): path to an expression matrix file
        network (str): path to a consensus network file
        out_dir (str): path to an output directory
        subnet (str, optional): path to a gene symbol file
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
              "spearman", "slope", "p-value")

    # Read in subnetwork file
    subnet_list = []
    out_subnet = None
    if subnet:
        out_subnet = open(pathlib.PurePath(out_dir).joinpath("sub_" + out_file_name + ".txt"), "w")
        with open(subnet) as subnet_file:
            for _id in subnet_file:
                subnet_list.append(_id.split("\n")[0].strip())

    exp = pd.read_csv(exp_mat, sep="\t", index_col=0)

    exp_dict = dict()
    with open(network, 'r') as fnet:
        fnet.readline()
        with open(pathlib.PurePath(out_dir).joinpath(out_file_name + ".txt"), 'w') as fout:
            fout.write('\t'.join(header) + '\n')
            for line in fnet:
                tokens = line.split('\t')
                node1 = tokens[0]
                node2 = tokens[1]
                mi = float(tokens[2])

                if node1 in exp_dict:
                    temp = exp_dict[node1]
                    gene_symbol1 = temp[0]
                    exp_values1 = temp[1]
                else:
                    exp_symbol_values1 = exp.loc[node1].values
                    gene_symbol1 = str(exp_symbol_values1[0])
                    exp_values1 = exp_symbol_values1[1:].astype(float)
                    exp_dict[node1] = (gene_symbol1, exp_values1)

                if node2 in exp_dict:
                    temp = exp_dict[node2]
                    gene_symbol2 = temp[0]
                    exp_values2 = temp[1]
                else:
                    exp_symbol_values2 = exp.loc[node2].values
                    gene_symbol2 = str(exp_symbol_values2[0])
                    exp_values2 = exp_symbol_values2[1:].astype(float)
                    exp_dict[node2] = (gene_symbol2, exp_values2)

                slope, intercept, r, p, stderr = stats.linregress(exp_values1, exp_values2)
                scc, sp = stats.spearmanr(exp_values1, exp_values2)
                pcc, pp = stats.pearsonr(exp_values1, exp_values2)

                row = (node1, node2, gene_symbol1, gene_symbol2,
                       "{0:.4f}".format(mi), "{0:.4f}".format(pcc), "{0:.4f}".format(scc),
                       "{0:.4f}".format(slope), "{0:.4f}".format(p))
                fout.write('\t'.join(row) + '\n')

                if out_subnet is not None:
                    if gene_symbol1 in subnet_list or gene_symbol2 in subnet_list:
                        out_subnet.write('\t'.join(row) + '\n')

    if out_subnet is not None:
        out_subnet.close()


def uprob(n):
    """ Implemented in statistics.py module inspired by Statistics::Distributions::uprob function in perl.
    Args:
        n (float): z-score
    Returns:
        p (float): p value
    """
    p = 0
    if abs(n) < 1.9:
        p = (1 + abs(n) * (0.049867347 + abs(n) * (
            0.0211410061 + abs(n) * 0.0032776263 + abs(n)
            * (0.0000380036 + abs(n) * (0.0000488906 + abs(n) * 0.000005383))))) ** (-16) / 2
    elif abs(n) <= 100:
        for i in range(18, 0, -1):
            p = i / (abs(n) + p)
        p = math.exp(-0.5 * abs(n) * abs(n)) / math.sqrt(2 * math.pi) / (abs(n) + p)
    if n < 0:
        p = 1 - p
    return p


if __name__ == '__main__':
    main()
