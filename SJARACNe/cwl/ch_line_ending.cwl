#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
doc: change Windows/Mac file line ending to Unix line ending

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_file)

baseCommand: ch_line_ending.py

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
      prefix: -i
      valueFrom: $(self.basename)
  output_file:
    type: string?
    inputBinding:
      position: 2
      prefix: -o

outputs:
  out_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_file == null) {
            return inputs.input_file.basename;
          }
          else {
            return inputs.output_file;
          }
        }
