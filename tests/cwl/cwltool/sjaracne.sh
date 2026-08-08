#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/sjaracne
cwltool --outdir ./results/cwl/sjaracne ./SJARACNe/cwl/sjaracne.cwl ./tests/inputs/cwl/brca_sjaracne.yml
