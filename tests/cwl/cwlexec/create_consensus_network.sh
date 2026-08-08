#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/cwlexec/create_consensus_network
cwlexec -pe PATH -c ./tests/cwl/cwlexec/config.json --outdir ./results/cwl/cwlexec/create_consensus_network \
./SJARACNe/cwl/create_consensus_network.cwl ./tests/inputs/cwl/brca_create_consensus_network.yml
