#!/usr/bin/env bash

cwltool --outdir ./test_data/outputs/cwl ./SJARACNe/cwl/create_consensus_network.cwl \
./test_data/inputs/cwl/create_consensus_network.yml
