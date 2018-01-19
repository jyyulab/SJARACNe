import sys


out_0 = open('00_cleanup_' + sys.argv[4] + '.sh', 'w')
out_1 = open('01_prepare_' + sys.argv[4] + '.sh', 'w')
out_2 = open('02_bootstrap_' + sys.argv[4] + '.sh', 'w')
out_3 = open('03_getconsensusnetwork_' + sys.argv[4] + '.sh', 'w')
thr = '1e-5'

exps =  open(sys.argv[1], 'r')
tfs = open(sys.argv[2], 'r')
sigs = open(sys.argv[3], 'r')

for tf in tfs:
        tf_dir = tf.strip().split('/')[0:-1]
        out_0.write('rm -rf ' + '/'.join(tf_dir) + '/sjaracne2_out_/ \n')
        out_0.write('rm -rf ' + '/'.join(tf_dir) + '/sjaracne2_out_.final \n')
        out_0.write('rm -rf ' + '/'.join(tf_dir) + '/sjaracne2_log_/ \n')
	out_1.write('perl -pe \'s/\\r\\n|\\n|\\r/\\n/g\' ' + tf.strip() + ' > ' + tf.strip() + '.tmp\n')
	out_1.write('rm ' + tf.strip() + '\n')
	out_1.write('mv ' + tf.strip() + '.tmp ' + tf.strip() + '\n') 

for sig in sigs:
        sig_dir = sig.strip().split('/')[0:-1]
        out_0.write('rm -rf ' + '/'.join(sig_dir) + '/sjaracne2_out_/ \n')
        out_0.write('rm -rf ' + '/'.join(sig_dir) + '/sjaracne2_out_.final \n')
        out_0.write('rm -rf ' + '/'.join(sig_dir) + '/sjaracne2_log_/ \n')
	out_1.write('perl -pe \'s/\\r\\n|\\n|\\r/\\n/g\' ' + sig.strip() + ' > ' + sig.strip() + '.tmp\n')
        out_1.write('rm ' + sig.strip() + '\n')
        out_1.write('mv ' + sig.strip() + '.tmp ' + sig.strip() + '\n')

for exp in exps:
	out_1.write('perl -pe \'s/\\r\\n|\\n|\\r/\\n/g\' ' + exp.strip() + ' > ' + exp.strip() + '.tmp\n')
        out_1.write('rm ' + exp.strip() + '\n')
        out_1.write('mv ' + exp.strip() + '.tmp ' + exp.strip() + '\n')

out_0.close()
out_1.close()

exps =  open(sys.argv[1], 'r')
tfs = open(sys.argv[2], 'r')
sigs = open(sys.argv[3], 'r')

for exp in exps:
	tf = tfs.readline().strip()
	if tf != "":
		tf_dir = tf.split('/')[0:-1]
		out_2.write('python bootstrap.py \"\" ' + '/'.join(tf_dir) + '/ 1 101 cluster \"-i ' + exp.strip() + ' -l ' + tf + ' -s ' + tf + ' -p 1e-7 -e 0 -a adaptive_partitioning -r 1 -H config/ -N ' + sys.argv[5] + '\"\nsleep 5\n')
		out_3.write('bsub -q compbio -P gn -M 16000 -o ' + '/'.join(tf_dir) + '/sjaracne2_log_/consensus_network.log python getconsensusnetwork.py ' + '/'.join(tf_dir) + '/sjaracne2_out_/ ' + thr + ' ' + '/'.join(tf_dir) + '/sjaracne2_out_.final/ ' + tf + '\n')
	sig = sigs.readline().strip()
	if sig != "":
		sig_dir = sig.split('/')[0:-1]
		out_2.write('python bootstrap.py \"\" ' + '/'.join(sig_dir) + '/ 1 101 cluster \"-i ' + exp.strip() + ' -l ' + sig + ' -s ' + sig + ' -p 1e-7 -e 0 -a adaptive_partitioning -r 1 -H config/ -N ' + sys.argv[5] + '\"\nsleep 5\n\n')
		out_3.write('bsub -q compbio -P gn -M 16000 -o ' + '/'.join(sig_dir) + '/sjaracne2_log_/consensus_network.log python getconsensusnetwork.py ' + '/'.join(sig_dir) + '/sjaracne2_out_/ ' + thr + ' ' + '/'.join(sig_dir) + '/sjaracne2_out_.final/ ' + sig + '\n\n')
		

out_2.close()
out_3.close()
