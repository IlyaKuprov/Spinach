# kernel/utilities/wigner_6j.m

- Signature: `w=wigner_6j(j1,j2,j3,j4,j5,j6)`

## Purpose

Wigner 6j-symbols. Syntax: w=wigner_6j(j1,j2,j3,j4,j5,j6) If physically inadmissible indices are supplied, a zero is returned. Order of elements: / j1 j2 j3 \ \ j4 j5 j6 /

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- j1-j6 -integers arranged in the order shown above

## Outputs

- w -the resulting 6j-symbol

## Implementation structure

- Wigner 6j-symbols. Syntax:
- w=wigner_6j(j1,j2,j3,j4,j5,j6)
- If physically inadmissible indices are supplied, a zero is
- returned. Order of elements:
- / j1 j2 j3 \
- \ j4 j5 j6 /
- j1-j6 -integers arranged in the order shown above
- w -the resulting 6j-symbol
- Check consistency
- Start from zero
- Use the definition
- Get the power index for -1
