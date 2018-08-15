sh data/sjaracne_BRCA100_scripts_/02_bootstrap_BRCA100.sh
jobs=$(ps -ef | grep "BRCA100" | grep "sjaracne -i" -c)
while [ $jobs -gt 0 ]
do
	sleep 30
	jobs=$(ps -ef | grep "BRCA100" | grep "sjaracne -i" -c)
done
sh data/sjaracne_BRCA100_scripts_/03_getconsensusnetwork_BRCA100.sh
jobs=$(ps -ef | grep "BRCA100" | grep getconsensusnetwork -c)
while [ $jobs -gt 0 ]
do
	sleep 30
	jobs=$(ps -ef | grep "BRCA100" | grep getconsensusnetwork -c)
done
sh data/sjaracne_BRCA100_scripts_/04_getenhancedconsensusnetwork_BRCA100.sh
