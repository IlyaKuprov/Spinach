# kernel/summaries/summary_couplings.m

- Signature: `summary_couplings(spin_system,header)`

## Purpose

Prints a summary of significant spin-spin coupling tensors in a Spinach spin system.

## Physical / mathematical content

For each reported pair, the table shows the coupling matrix components and its isotropic contribution, calculated as the matrix trace divided by three. It also reports the norms of the rank-1 and rank-2 parts obtained from the spherical-tensor decomposition.

## Numerical / algorithmic content

The routine selects entries whose matrix 2-norm exceeds `2*pi*spin_system.tols.inter_cutoff`. Matrix components, the isotropic contribution, and rank norms are reported after division by `2*pi`; the rank norms use the Frobenius norm. Input checks require a structure and a character-string header.

## Parameters / inputs

- spin_system - Spinach spin system description object
- header - a string of text to precede the summary

## Outputs

- Prints the summary through `report.m` to the console or user-specified output.

## Implementation structure

After printing column headings, the function traverses the selected matrix entries and reports the spin indices, matrix elements, isotropic value, and rank-1 and rank-2 norms. It ends the table with a separator.
