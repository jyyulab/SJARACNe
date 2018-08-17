import os, sys

SJARACNE_PATH = os.environ['SJARACNE_PATH']
SJARACNE_PATH = SJARACNE_PATH if SJARACNE_PATH.endswith('/') else SJARACNE_PATH + '/'
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, SJARACNE_PATH)
sys.path.insert(0, PYTHON_PATH)

if sys.platform == 'linux':
	sjaracne = 'sjaracne'
elif sys.platform == 'darwin':
	sjaracne = 'sjaracne.osx'

import argparse
import numpy as np
import shlex, subprocess

def setup(args):
	parser = argparse.ArgumentParser(description='SJARACNE package.')
	parser.add_argument('project_name', help='Project name')
	parser.add_argument('expression_matrix', help='Path to expression matrix, row indexes are used as the nodes in the network.')
	parser.add_argument('hub_genes', help='Path to hub genes, containing list of symbols to be considered as hub genes.')
	parser.add_argument('--bootstrap', type=int, default=100, help='Number of bootstrap networks.')
	parser.add_argument('--c_threshold', type=float, default=float('1e-5'), help='P-value threshold in building consensus network.')
	parser.add_argument('--p_threshold', type=float, default=float('1e-7'), help='P-value threshold in building bootstrap netwroks.')
	parser.add_argument('--depth', type=int, default=40, help='Maximum partitioning depth.')
	parser.add_argument('--run', type=bool, default=False, help='Whether run the pipeline or just generate and stop.')
	parser.add_argument('--host', default='LOCAL', help='Whether to run on clusters or localhost. [LOCAL | CLUSTER]')
	parser.add_argument('outdir', help='Output directory')
	args_ = parser.parse_args(args[1:])
	if args_.host == 'CLUSTER':
		args_.run = 'False'
	return args_

def setup_directory(out_dir, project_name):
	boot_path = out_dir + 'sjaracne_' + project_name + '_out_/'
	log_path = out_dir + 'sjaracne_' + project_name + '_log_/'
	net_path = out_dir + 'sjaracne_' + project_name + '_out_.final/'
	script_path = out_dir + 'sjaracne_' + project_name + '_scripts_/'
	if not os.path.exists(boot_path):
		os.mkdir(boot_path)
	if not os.path.exists(log_path):
		os.mkdir(log_path)
	if not os.path.exists(net_path):
		os.mkdir(net_path)
	if not os.path.exists(script_path):
		os.mkdir(script_path)
	return [boot_path, log_path, net_path, script_path]

def cleanup(args, paths):
	out_0 = open(paths[3] + '00_cleanup_' + args.project_name + '.sh', 'w')
	out_0.write('rm -rf ' + paths[0] + '\n')
	out_0.write('rm -rf ' + paths[1] + '\n')
	out_0.write('rm -rf ' + paths[2] + '\n')
	out_0.write('rm -rf ' + paths[3] + '\n')
	out_0.close()

def prep(args, paths):
	out_1 = open(paths[3] + '01_prepare_' + args.project_name + '.sh', 'w')
	out_1.write('perl -pe \'s/\\r\\n|\\n|\\r/\\n/g\' ' + args.hub_genes + ' > ' + args.hub_genes + '.tmp\n')
	out_1.write('rm ' + args.hub_genes + '\n')
	out_1.write('mv ' + args.hub_genes + '.tmp ' + args.hub_genes + '\n')
	out_1.write('perl -pe \'s/\\r\\n|\\n|\\r/\\n/g\' ' + args.expression_matrix + ' > ' + args.expression_matrix + '.tmp\n')
	out_1.write('rm ' + args.expression_matrix + '\n')
	out_1.write('mv ' + args.expression_matrix + '.tmp ' + args.expression_matrix + '\n')
	out_1.close()

def bootstrap(args, paths):
	b = args.bootstrap + 1
	out_2 = open(paths[3] + '02_bootstrap_' + args.project_name + '.sh', 'w')
	for i in np.arange(1, b):
		fname = paths[0] + args.project_name + '_run_' + str(i).zfill(int(np.log10(b)) + 1) + '.adj'
		lname = paths[1] + args.project_name + '_run_' + str(i).zfill(int(np.log10(b)) + 1) + '.log'
		arg = ' -i ' + args.expression_matrix + ' -l ' + args.hub_genes + ' -s ' + args.hub_genes + ' -p ' + str(args.p_threshold) + ' -e 0 -a adaptive_partitioning -r 1 -H ' + SJARACNE_PATH + 'config/ -N ' + str(args.depth)
		script = SJARACNE_PATH + 'bin/' + sjaracne + ' ' + arg + ' -o ' + fname + ' -S ' + str(i) + ' '
		if args.host == 'LOCAL':
			script += ' >> ' + lname + ' &'
		script += '\n'
		out_2.write(script)
	out_2.close()

def consensus(args, paths):
	out_3 = open(paths[3] + '03_getconsensusnetwork_' + args.project_name + '.sh', 'w')
	script = PYTHON_PATH + ' ' + SJARACNE_PATH + 'getconsensusnetwork.py ' + paths[0] + ' ' + str(args.c_threshold) + ' ' + paths[2] + ' '
	if args.host == 'LOCAL':
		script += ' >> ' + paths[1] + args.project_name + '_consensus_network.log'
	script += '\n'
	out_3.write(script)
	out_3.close()

def enhanced(args, paths):
	out_4 = open(paths[3] + '04_getenhancedconsensusnetwork_' + args.project_name + '.sh', 'w')
	script = PYTHON_PATH + ' ' + SJARACNE_PATH + 'getenhancedconsensusnetwork.py ' + args.expression_matrix + ' ' + paths[2] + 'consensus_network_3col_.txt ' + paths[2] + ' '
	if args.host == 'LOCAL':
		script += ' >> ' + paths[1] + args.project_name + '_enhanced_network.log'
	script += '\n'
	out_4.write(script)
	out_4.close()

def pipeline(args, paths):
	out_0 = open(paths[3] + '00_pipeline_' + args.project_name + '.sh', 'w')
	#script = 'sh ' + paths[3] + '00_cleanup_' + args.project_name + '.sh\n'
	#out_0.write(script)
	#script = 'sh ' + paths[3] + '01_prepare_' + args.project_name + '.sh\n'
	#out_0.write(script)
	script = 'sh ' + paths[3] + '02_bootstrap_' + args.project_name + '.sh\n'
	out_0.write(script)
	out_0.write('jobs=$(ps -ef | grep \"' + args.project_name + '\" | grep \"' + sjaracne + ' -i\" -c)\n')
	out_0.write('while [ $jobs -gt 0 ]\ndo\n\tsleep 30\n\tjobs=$(ps -ef | grep \"' + args.project_name + '\" | grep \"' + sjaracne + ' -i\" -c)\ndone\n')
	script = 'sh ' + paths[3] + '03_getconsensusnetwork_' + args.project_name + '.sh\n'
	out_0.write(script)
	out_0.write('jobs=$(ps -ef | grep \"' + args.project_name + '\" | grep getconsensusnetwork -c)\n')
	out_0.write('while [ $jobs -gt 0 ]\ndo\n\tsleep 30\n\tjobs=$(ps -ef | grep \"' + args.project_name + '\" | grep getconsensusnetwork -c)\ndone\n')
	script = 'sh ' + paths[3] + '04_getenhancedconsensusnetwork_' + args.project_name + '.sh\n'
	out_0.write(script)
	out_0.close()

def run(args):
	args_ = setup(args)
	paths = setup_directory(args_.outdir, args_.project_name)
	cleanup(args_, paths)
	prep(args_, paths)
	bootstrap(args_, paths)
	consensus(args_, paths)
	enhanced(args_, paths)
	pipeline(args_, paths)
	if args_.run == True:
		script = 'sh ' + paths[3] + '00_pipeline_' + args_.project_name + '.sh '
		if args_.host == 'LOCAL':
			script += ' >> ' + paths[1] + args_.project_name + '_pipeline.log'
		script += '\n'
		subprocess.Popen(shlex.split(script))

if __name__ == '__main__':
	run(sys.argv)
	print('[INFO] --> [ARCN] Finished.')
