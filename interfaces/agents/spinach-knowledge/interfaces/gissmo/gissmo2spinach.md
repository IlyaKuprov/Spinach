# interfaces/gissmo/gissmo2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/gissmo/gissmo2spinach.m`
- Signature: `[sys,inter]=gissmo2spinach(filename,subsystem)`
- Total lines: 206

## Purpose

Reads GISSMO files and forms Spinach data structures. Syntax: [sys,inter]=gissmo2spinach(file_name,subsystem)

## Physical / mathematical content

- This file belongs to the `interfaces` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- file_name -character string with the name of
- the GISSMO XML file
- subsystem -which of the coupling matrices to
- to import

## Outputs

- sys, inter -Spinach data structures, ready for
- calling create.m
- Note: GISSMO only provides chemical shifts, J-couplings, the
- non-selective line width, and the magnet field. You may
- want to add further parameters by editing sys and inter
- data structures manually.

## Implementation structure

- Reads GISSMO files and forms Spinach data structures. Syntax:
- [sys,inter]=gissmo2spinach(file_name,subsystem)
- file_name - character string with the name of
- the GISSMO XML file
- subsystem - which of the coupling matrices to
- to import
- sys, inter - Spinach data structures, ready for
- calling create.m
- Note: GISSMO only provides chemical shifts, J-couplings, the
- non-selective line width, and the magnet field. You may
- want to add further parameters by editing sys and inter
- data structures manually.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `parsexml()`, `false()`, `strcmpi()`, `str2double()`, `spin()`, `true()`, `fwhm2rlx()`, `ischar()`, `exist()`.
