psub -K -P TF -J TF_SJARACNE_Bootstrap -q compbio -M 2000 -i ./test_data/outputs/SJARACNE_TF/SJARACNE_scripts/02_bootstrap_TF.sh -oo ./test_data/outputs/SJARACNE_TF/SJARACNE_log/TF_SJARACNE_Bootstrap.%J.%I.out -eo ./test_data/outputs/SJARACNE_TF/SJARACNE_log/TF_SJARACNE_Bootstrap.%J.%I.err 
sleep 30
psub -K -P TF -J TF_SJARACNE_Consensus -q compbio -M 2000 -i ./test_data/outputs/SJARACNE_TF/SJARACNE_scripts/03_getconsensusnetwork_TF.sh -oo ./test_data/outputs/SJARACNE_TF/SJARACNE_log/TF_SJARACNE_Consensus.%J.%I.out -eo ./test_data/outputs/SJARACNE_TF/SJARACNE_log/TF_SJARACNE_Consensus.%J.%I.err 
sleep 30
psub -K -P TF -J TF_SJARACNE_Enhanced -q compbio -M 2000 -i ./test_data/outputs/SJARACNE_TF/SJARACNE_scripts/04_getenhancedconsensusnetwork_TF.sh -oo ./test_data/outputs/SJARACNE_TF/SJARACNE_log/TF_SJARACNE_Enhanced.%J.%I.out -eo ./test_data/outputs/SJARACNE_TF/SJARACNE_log/TF_SJARACNE_Enhanced.%J.%I.err 
