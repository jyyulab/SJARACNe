#!/usr/bin/env bash

export PATH="$(pwd)/SJARACNe/bin:$PATH"
mkdir -p ./results/cwl/create_consensus_network
cwltool --outdir ./results/cwl/create_consensus_network ./SJARACNe/cwl/create_consensus_network.cwl \
./tests/inputs/cwl/brca_create_consensus_network.yml
