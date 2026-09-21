# interfaces/spinxml/x2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/spinxml/x2spinach.m`
- Signature: `[sys,inter]=x2spinach(filename,shielding_refs)`
- Total lines: 895

## Purpose

Reads SpinXML files and forms Spinach data structures. Syntax: [sys,inter]=x2spinach(file_name,shielding_refs)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `parsexml()`, `strcmpi()`, `str2double()`, `xyz()`, `clear()`, `isequal()`, `spin_ids()`, `anas2mat()`, `axrh2mat()`, `spsk2mat()`, `euler2dcm()`, `qter2dcm()`, `dcm()`, `anax2dcm()`, `iselectron()`.
