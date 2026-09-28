# kernel/summaries/summary_zeeman.m

- Signature: `summary_zeeman(spin_system,header)`

## Purpose

Prints a summary of the Zeeman interaction matrices stored in a Spinach system.

## Parameters / inputs

- `spin_system` - Spinach spin system structure.
- `header` - character string printed before the table.

## Output

No MATLAB output argument; writes through `report`.

## Behavior

The function checks the inputs, then processes each nonempty Zeeman matrix. For each, it reports the spin index, isotope, multiplicity (2S+1), matrix, isotropic part (trace/3), and the matrix 2-norms of the rank-1 and rank-2 components. It obtains those components with `mat2sphten` and `sphten2mat`, and prints the third row of the original matrix as well.
