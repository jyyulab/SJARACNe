#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/cwlexec/ch_line_ending
cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./results/cwl/cwlexec/ch_line_ending ./SJARACNe/cwl/ch_line_ending.cwl ./tests/inputs/cwl/brca_ch_line_ending.yml
