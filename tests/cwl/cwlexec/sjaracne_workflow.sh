#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/cwlexec/sjaracne_workflow
cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./results/cwl/cwlexec/sjaracne_workflow ./SJARACNe/cwl/sjaracne_workflow.cwl ./tests/inputs/cwl/brca_sjaracne_workflow.yml
