#!/usr/bin/env bash

cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./test_data/outputs/cwl \
./SJARACNe/cwl/create_consensus_network.cwl ./test_data/inputs/cwl/create_consensus_network.yml
