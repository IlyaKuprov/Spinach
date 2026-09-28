# kernel/summaries/summary_basis.m

- Signature: `summary_basis(spin_system)`

## Purpose

Prints a summary of the basis set for a Spinach system. Syntax: `summary_basis(spin_system)`.

## Physical / mathematical content

For each reported basis state, lists the irreducible spherical-tensor quantum-number pairs `(L,M)` for each spin.

## Numerical / algorithmic content

Reports the basis dimension and its percentage of the full state space. If the number of basis states exceeds `spin_system.tols.basis_hush`, detailed state labels are suppressed.

## Parameters / inputs

- `spin_system` - Spinach spin system description object.

## Outputs

- Prints through `report.m` to the console or the user-specified output.

## Implementation structure

- Checks that the input is a structure, gets the basis dimension, conditionally reports each basis state's `(L,M)` labels, then reports the dimension and percentage of the full state space.
