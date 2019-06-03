#!/usr/bin/env bash

export PATH=`pwd`/SJARACNe/bin:$PATH
./SJARACNe/sjaracne.py local \
-e ./test_data/inputs/BRCA100.exp \
-g ./test_data/inputs/tf.txt \
-n 2 \
-o ./test_data/outputs/cwl/cwltool/SJARACNE_out.final
