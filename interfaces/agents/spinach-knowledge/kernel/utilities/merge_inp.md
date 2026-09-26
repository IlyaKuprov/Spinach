# kernel/utilities/merge_inp.m

- Signature: `[sys,inter]=merge_inp(sys_parts,inter_parts)`

## Purpose

Merges multiple sys and inter structures into one. Useful for setting up chemical kinetics simulations where the molecules come from different DFT calculations. Syntax: [sys,inter]=merge_inp(sys_parts,inter_parts)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

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

- Merges multiple sys and inter structures into one. Useful for
- setting up chemical kinetics simulations where the molecules
- come from different DFT calculations. Syntax:
- [sys,inter]=merge_inp(sys_parts,inter_parts)
- sys_parts -a cell array of sys structures
- to be merged
- inter_parts -a cell array of inter structures
- sys -resulting sys structure
- inter -resulting inter structure
- Note: extensive fields are concatenated with spin and chemical
- subsystem indices offset as appropriate; non-extensive
- fields (magnet, temperature, relaxation settings, etc.)
