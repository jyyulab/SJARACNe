#!/usr/bin/env bash

export PATH=`pwd`/SJARACNe/bin:$PATH
cwltool --outdir ./test_data/outputs/cwl ./SJARACNe/cwl/ch_line_ending.cwl ./test_data/inputs/cwl/ch_line_ending.yml