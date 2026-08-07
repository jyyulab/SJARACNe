#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/cwlexec/sjaracne
cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./results/cwl/cwlexec/sjaracne ./SJARACNe/cwl/sjaracne.cwl ./tests/inputs/cwl/brca_sjaracne.yml
