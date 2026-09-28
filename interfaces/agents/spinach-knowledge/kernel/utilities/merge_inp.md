# kernel/utilities/merge_inp.m

- Signature: `[sys,inter]=merge_inp(sys_parts,inter_parts)`

## Purpose

Merges multiple sys and inter structures into one. Useful for setting up chemical kinetics simulations where the molecules come from different DFT calculations. Syntax: [sys,inter]=merge_inp(sys_parts,inter_parts)

## Physical / mathematical content

- Combines `sys` and `inter` structures for multiple subsystems, including supported interaction groups and chemical-process specifications.

## Numerical / algorithmic content

- Concatenates extensive subsystem data, offsets spin-index lists by preceding subsystem spin counts, and combines supported square coupling arrays block-diagonally. Shared values must agree across subsystems; inconsistent or unhandled fields cause errors.

## Parameters / inputs

- sys_parts -a cell array of sys structures
- to be merged
- inter_parts -a cell array of inter structures
- to be merged

## Outputs

- sys -resulting sys structure
- inter -resulting inter structure
- Note: extensive fields are concatenated with spin and chemical
- subsystem indices offset as appropriate; non-extensive
- fields (magnet, temperature, relaxation settings, etc.)
- must be identical in all subsystems that supply them, and
- any difference is treated as an error. Nested groups
- (zeeman, coupling, giant, suscept, chem) and every field
- must be present in all subsystems or in none. An error is
- thrown for unhandled subfields. Coordinates and suscepti-
- bility centres from all subsystems are assumed to refer
- to one common frame of reference; spin index lists are
- returned as row vectors.

## Implementation structure

- Validates paired, non-empty row cell arrays and counts the spins in each subsystem.
- Merges common fields, per-spin arrays, coordinates, coupling blocks, spin-index lists, and supported nested groups with field-specific handlers.
- Rejects partial nested groups, unequal shared values, and any unhandled fields.
