# kernel/utilities/isoswap.m

- Signature: `[sys,inter]=isoswap(sys,inter,spins,new_iso)`

## Purpose

Makes isotope replacements in the input structures. All interactions are automatically scaled as necessary. Syntax: [sys,inter]=isoswap(sys,inter,spins,new_iso)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- sys, inter -Spinach input data structures
- spins -a vector of integers specifying
- spin numbers
- new_iso -character string specifying
- the new isotope
- Output:
- sys, inter -Spinach input data structures
- Note: quadratic and higher order couplings are wiped and a warning
- is printed -those are not transferable.

## Implementation structure

- Makes isotope replacements in the input structures. All interactions
- are automatically scaled as necessary. Syntax:
- [sys,inter]=isoswap(sys,inter,spins,new_iso)
- sys, inter -Spinach input data structures
- spins -a vector of integers specifying
- spin numbers
- new_iso -character string specifying
- the new isotope
- Output:
- Note: quadratic and higher order couplings are wiped and a warning
- is printed -those are not transferable.
- Grumbler missing
