#!/usr/bin/env bash

export PATH=`pwd`/SJARACNe/bin:$PATH
cwltool --parallel --outdir ./test_data/outputs/cwl/cwltool ./SJARACNe/cwl/sjaracne_workflow.cwl ./test_data/inputs/cwl/sjaracne_workflow.yml
