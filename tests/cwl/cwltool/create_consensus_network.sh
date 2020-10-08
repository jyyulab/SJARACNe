#!/usr/bin/env bash

export PATH=`pwd`/SJARACNe/bin:$PATH
cwltool --outdir ./test_data/outputs/cwl ./SJARACNe/cwl/create_consensus_network.cwl \
./test_data/inputs/cwl/create_consensus_network.yml
