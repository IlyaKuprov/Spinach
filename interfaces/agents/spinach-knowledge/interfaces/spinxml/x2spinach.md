# interfaces/spinxml/x2spinach.m

- Signature: `[sys,inter]=x2spinach(filename,shielding_refs)`

## Purpose

Reads SpinXML files and forms Spinach data structures. Syntax: [sys,inter]=x2spinach(file_name,shielding_refs)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- shielding_refs -when absolute shielding is provided, this
- array must give absolute shielding of the
- corresponding nuclei in the reference sub-
- stance, e.g. {{'1H',31.5},{'13C',189.7}}.
- This is necessary because Spinach requires
- chemical shifts rather than shieldings.
- If chemical shifts are provided, use an
- empty cell array.

## Outputs

- sys, inter -Spinach spin system input structures
- WARNING: this function assumes that the SpinXML file has passed
- the validation against the schema, which may be obtain-
- ed from http://spindynamics.org/SpinXML.php

## Implementation structure

- Reads SpinXML files and forms Spinach data structures. Syntax:
- [sys,inter]=x2spinach(file_name,shielding_refs)
- shielding_refs -when absolute shielding is provided, this
- array must give absolute shielding of the
- corresponding nuclei in the reference sub-
- stance, e.g. {{'1H',31.5},{'13C',189.7}}.
- This is necessary because Spinach requires
- chemical shifts rather than shieldings.
- If chemical shifts are provided, use an
- empty cell array.
- sys, inter -Spinach spin system input structures
- WARNING: this function assumes that the SpinXML file has passed
