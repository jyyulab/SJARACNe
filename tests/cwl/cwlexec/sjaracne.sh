#!/usr/bin/env bash

export PATH=`pwd`/SJARACNe/bin:$PATH
cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./test_data/outputs/cwl ./SJARACNe/cwl/sjaracne.cwl ./test_data/inputs/cwl/sjaracne.yml
