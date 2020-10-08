#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: validation of input files

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.exp_file)
      - $(inputs.probe_file)

baseCommand: QC_input.py

inputs:
  exp_file:
    type: File
    inputBinding:
      position: 1
      prefix: -e
      valueFrom: $(self.basename)
  probe_file:
    type: File
    inputBinding:
      position: 2
      prefix: -g
      valueFrom: $(self.basename)

outputs: []
