#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: get_consensus_network.py

inputs:
  adjmat_dir:
    type: Directory
    inputBinding:
      position: 1
    doc: directory with adjacent matrix

  p_thresh_arg:
    type: float
    inputBinding:
      position: 2
    doc: P value threshold

  output_dir:
    type: string
    inputBinding:
      position: 3
    doc: output directory name

outputs:
  out_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)
