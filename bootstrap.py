##############################################################
## Python script to bootstrap multiple jobs and submit them ##
## to computing clusters at St. Jude Children's Research    ##
## Hospital.                                                ##
## Based on the given arguments, this script creates        ##
## multiple shell script and submit each shell script to    ##
## the computing clusters.                                  ##
## Author: Alireza Khatamian				   				##
## Email: akhatami@stjude.org				    			##
##############################################################
import os, sys, shlex, subprocess

if len(sys.argv) < 5:
	print "Insufficient parameters!"
	print "Usage: python bootstrap.py <output_prefix> <out_dest> <lower#> <upper#> <local | cluster> <sjaracne2 arguments>"
	print "Example: python bootstrap.py out_prefix input/ 1 101 local \"-i input/data.exp -l input/tflist.txt -s input/tflist.txt -p 1e-7 -e 0 -a adaptive_partitioning -r 1\""
	exit()
# Creating and submitting shell scripts for each run with respect to lower and upper bound arguments
os.mkdir(sys.argv[2] + 'sjaracne2_out_' + sys.argv[1])
os.mkdir(sys.argv[2] + 'sjaracne2_out_' + sys.argv[1] + '.final')
os.mkdir(sys.argv[2] + 'sjaracne2_log_' + sys.argv[1])
for i in range(int(sys.argv[3]), int(sys.argv[4])):
	fname = sys.argv[2] + 'sjaracne2_out_' + sys.argv[1] + '/' + sys.argv[1] + '_run_' + str(i).zfill(3) + '.adj'
	lname = sys.argv[2] + 'sjaracne2_log_' + sys.argv[1] + '/' + sys.argv[1] + '_run_' + str(i).zfill(3) + '.log'
	if sys.argv[5] == 'local':
		script = './sjaracne ' + sys.argv[6] + ' -o ' + fname + ' -S ' + str(i) + ' > ' + lname
	elif sys.argv[5] == 'cluster':
		script = 'bsub -q compbio -P gn -o ' + lname + ' -R \"rusage[mem=16000]\" ./sjaracne ' + sys.argv[6] + ' -o ' + fname + ' -S ' + str(i)
	# Submit the job within the shell script into the computing clusters
	subprocess.Popen(shlex.split(script))
