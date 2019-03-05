sh ./test_data/outputs/local/SJARACNE_TF/SJARACNE_scripts/02_bootstrap_TF.sh
jobs=$(ps -ef | grep "TF" | grep sjaracne -c)
while [ $jobs -gt 0 ]
do
	sleep 30
	jobs=$(ps -ef | grep "TF" | grep sjaracne -c)
done
sh ./test_data/outputs/local/SJARACNE_TF/SJARACNE_scripts/03_getconsensusnetwork_TF.sh >> ./test_data/outputs/local/SJARACNE_TF/SJARACNE_log/TF_consensus_network.out
jobs=$(ps -ef | grep "TF" | grep getconsensusnetwork -c)
while [ $jobs -gt 0 ]
do
	sleep 30
	jobs=$(ps -ef | grep "TF" | grep getconsensusnetwork -c)
done
sh ./test_data/outputs/local/SJARACNE_TF/SJARACNE_scripts/04_getenhancedconsensusnetwork_TF.sh >> ./test_data/outputs/local/SJARACNE_TF/SJARACNE_log/TF_enhanced_network.out
