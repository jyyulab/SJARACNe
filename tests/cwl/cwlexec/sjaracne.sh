#!/usr/bin/env bash

cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./test_data/outputs/cwl ./SJARACNe/cwl/sjaracne.cwl ./test_data/inputs/cwl/sjaracne.yml
