#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/sjaracne_workflow
cwltool --parallel --outdir ./results/cwl/sjaracne_workflow ./SJARACNe/cwl/sjaracne_workflow.cwl ./tests/inputs/cwl/brca_sjaracne_workflow.yml
