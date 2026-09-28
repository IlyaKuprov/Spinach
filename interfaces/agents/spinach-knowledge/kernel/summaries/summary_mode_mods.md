# kernel/summaries/summary_mode_mods.m

- Signature: `summary_mode_mods(spin_system,header)`

## Purpose

Summarizes how dimensionless bosonic mode displacement coordinates modulate spin-spin coupling tensors and effective local fields in the spin Hamiltonian.

## Physical / mathematical content

For each reported modulation, the table identifies the mode pair and derivative order, names the affected coupling tensor or local field, and gives the derivative norm in Hz.

## Numerical / algorithmic content

The function traverses nonempty entries in `spin_system.inter.modes.coupling_mod` and `zeeman_mod`. Coupling-tensor derivatives use the Frobenius norm; effective-field derivatives use the 2-norm. Both are divided by `2*pi` for reporting. Input checks require a structure and a character-string header.

## Parameters / inputs

- spin_system - Spinach spin system description object
- header - a string of text to precede the summary

## Outputs

- Prints the modulation summary through `report.m` to the console or user-specified output.

## Implementation structure

After validation and table headings, the function reports nonempty coupling-tensor derivatives and effective-field derivatives, with their mode indices, order, target, and norm.
