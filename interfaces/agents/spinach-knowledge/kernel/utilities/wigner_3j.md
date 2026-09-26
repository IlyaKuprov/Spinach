# kernel/utilities/wigner_3j.m

- Signature: `w=wigner_3j(j1,m1,j2,m2,j3,m3)`

## Purpose

Calculates Wigner 3j-symbols. Syntax: w=wigner_3j(j1,m1,j2,m2,j3,m3) If physically inadmissible indices are supplied, a zero is returned. Order of elements: /j1 j2 j3\ \m1 m2 m3/

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- j1-j3 -integers arranged in the order shown above
- m1-m3 -integers arranged in the order shown above

## Outputs

- w -the resulting 3j-symbol

## Implementation structure

- Calculates Wigner 3j-symbols. Syntax:
- w=wigner_3j(j1,m1,j2,m2,j3,m3)
- If physically inadmissible indices are supplied, a zero is
- returned. Order of elements:
- /j1 j2 j3\
- \m1 m2 m3/
- j1-j3 -integers arranged in the order shown above
- m1-m3 -integers arranged in the order shown above
- w -the resulting 3j-symbol
- Check consistency
- Call Clebsch-Gordan coefficients
- Consistency enforcement
