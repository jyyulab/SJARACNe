#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/ch_line_ending
cwltool --outdir ./results/cwl/ch_line_ending ./SJARACNe/cwl/ch_line_ending.cwl ./tests/inputs/cwl/brca_ch_line_ending.yml
