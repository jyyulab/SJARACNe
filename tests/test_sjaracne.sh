#!/usr/bin/env bash

./SJARACNe/sjaracne.py lsf \
-e /research/rgs01/home/clusterHome/lding/Git/SJARACNe/test_data/inputs/BRCA100.exp \
-g /research/rgs01/home/clusterHome/lding/Git/SJARACNe/test_data/inputs/tf.txt \
-c /home/lding/Git/SJARACNe/SJARACNe/config/ \
-n 2 \
-o ./test_data/outputs/cwl/cwltool/SJARACNE_out.final \
-q priority
