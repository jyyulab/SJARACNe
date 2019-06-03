#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: ExpressionTool

requirements:
 - class: InlineJavascriptRequirement

inputs:
  number:
    type: int
    label: a positive integer

outputs:
  str_array:
    type: string[]

expression: |
  ${ var str_arr = [], str_i = '';
     for (var i = 0; i < inputs.number; i++) {
       if (i <= 999) {
         str_i = ("00" + i).slice(-3);
       }
       str_arr.push('TF_run_' + str_i + '.adj');
     }
     return { "str_array": str_arr };
   }