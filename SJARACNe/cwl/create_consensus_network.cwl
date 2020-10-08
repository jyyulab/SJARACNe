#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: create_consensus_network.py

inputs:
  adjmat_dir:
    type: Directory
    inputBinding:
      position: 1
      prefix: -a
    doc: directory with adjacent matrix

  p_thresh_arg:
    type: float
    inputBinding:
      position: 2
      prefix: -p
    doc: P value threshold

  exp_mat:
    type: File
    inputBinding:
      position: 3
      prefix: -e
    doc: expression matrix file

  output_dir:
    type: string
    inputBinding:
      position: 4
      prefix: -o
    doc: output directory name

  subnet:
    type: File?
    inputBinding:
      position: 5
      prefix: -s
    doc: file with gene symbols of interest to build a subnet

outputs:
  out_dir:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/consensus_network_ncol_.txt
