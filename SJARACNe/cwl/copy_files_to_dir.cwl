#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: ExpressionTool

requirements:
 - class: InlineJavascriptRequirement

inputs:
  input_files:
    type: File[]
  dirname:
    type: string

outputs:
  out_dir:
    type: Directory

expression: |
  ${
    return {"out_dir": { "class": "Directory", "basename": inputs.dirname, "listing": inputs.input_files} };
  }
