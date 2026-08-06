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
  output_file:
    type: string
    default: hub_overlap_validation.txt
    inputBinding:
      position: 3
      prefix: -o

outputs:
  validation_report:
    type: File
    outputBinding:
      glob: $(inputs.output_file)
