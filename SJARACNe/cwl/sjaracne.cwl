#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: Scalable solution of ARACNe that dramatically improves the computational performance.

baseCommand: sjaracne.exe

inputs:
  exp_file:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    doc: Input gene expression profile dataset (required)
  probe_file_tf:
    type: File
    inputBinding:
      position: 2
      prefix: -l
    doc: File containing a list of probes annotated as transcription factors in the input dataset
  probe_file_subnetwork:
    type: File
    inputBinding:
      position: 3
      prefix: -s
    doc: File containing a list of probes for which a subnetwork will be constructed
  p_value:
    type: float
    inputBinding:
      position: 4
      prefix: -p
  tolerance:
    type: int
    default: 0
    inputBinding:
      position: 5
      prefix: -e
    doc: DPI tolerance
  algorithm:
    type: string
    default: adaptive_partitioning
    inputBinding:
      position: 6
      prefix: -a
    doc: algorithm
  sample_number:
    type: int
    default: 1
    inputBinding:
      position: 7
      prefix: -r
    doc: Bootstrap sample number
  aracne_config_dir:
    type: Directory
    inputBinding:
      position: 8
      prefix: -H
    doc: Directory containing ARACNe configuration files, default is current working directory
  npar_limit:
    type: int
    default: 20
    inputBinding:
      position: 9
      prefix: -N
    doc: Maximum allowed value of npar
  output_file_name:
    type: string
    inputBinding:
      position: 10
      prefix: -o
    doc: Output file name (optional)
  seed:
    type: int
    default: 1
    inputBinding:
      position: 11
      prefix: -S
    doc: Initial seed for random number generator

outputs:
  out_adj:
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)
