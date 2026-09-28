# kernel/summaries/summary_symmetry.m

- Signature: `summary_symmetry(spin_system,header)`

## Purpose

Prints the permutation-symmetry groups and their associated spins from a Spinach system.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.
- `header` - character string printed before the table.

## Output

No MATLAB output argument; writes the header and table through `report`.

## Behavior

The function checks that `spin_system` is a structure and `header` is a character string. It iterates over `spin_system.comp.sym_spins` and prints each corresponding entry from `spin_system.comp.sym_group` alongside its spin list.
