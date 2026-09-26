# interfaces/gissmo/gissmo2spinach.m

- Signature: `[sys,inter]=gissmo2spinach(filename,subsystem)`

## Purpose

Reads GISSMO files and forms Spinach data structures. Syntax: [sys,inter]=gissmo2spinach(file_name,subsystem)

## Physical / mathematical content

## Numerical / algorithmic content

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
