# kernel/utilities/gtensorof.m

- Signature: `g=gtensorof(spin_system,spin_number)`

## Purpose

Returns the g-tensor of the specified spin at the input orientation. Syntax: g=gtensorof(spin_system,spin_number)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- spin_number -a positive integer specifying the number
- of the spin in the sys.isotopes list

## Outputs

- g -a 3x3 matrix in Bohr magneton units
- Note: the same convention (mu=-mu_b*g*S/hbar) is used for
- the nuclei, meaning that their g-tensors are much
- smaller than those of electrons.

## Implementation structure

- Returns the g-tensor of the specified spin at the input
- orientation. Syntax:
- g=gtensorof(spin_system,spin_number)
- spin_number -a positive integer specifying the number
- of the spin in the sys.isotopes list
- g -a 3x3 matrix in Bohr magneton units
- Note: the same convention (mu=-mu_b*g*S/hbar) is used for
- the nuclei, meaning that their g-tensors are much
- smaller than those of electrons.
- Check consistency
- Compute the g-tensor
- Consistency enforcement
