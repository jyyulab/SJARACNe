#!/usr/bin/env bash

cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./test_data/outputs/cwl/cwlexec ./SJARACNe/cwl/sjaracne_workflow.cwl ./test_data/inputs/cwl/sjaracne_workflow.yml