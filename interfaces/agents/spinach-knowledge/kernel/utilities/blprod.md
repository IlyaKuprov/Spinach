# kernel/utilities/blprod.m

- Signature: `[X1_AB,X2_AB]=blprod(A,B)`

## Purpose

Extension of Blicharski's tensor invariants into scalar products of different spin interaction tensors using polarisation identi- ties. Syntax: [X1_AB,X2_AB]=blprod(A,B)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- A -a real 3x3 matrix
- B -a real 3x3 matrix

## Outputs

- X1_AB -cross-correlation amplitude, first rank
- X2_AB -cross-correlation amplitude, second rank
- Note: this function is not sensitive to the isotropic components
- of A and B tensors.

## Implementation structure

- Extension of Blicharski's tensor invariants into scalar products
- of different spin interaction tensors using polarisation identi-
- ties. Syntax:
- [X1_AB,X2_AB]=blprod(A,B)
- A - a real 3x3 matrix
- B - a real 3x3 matrix
- X1_AB - cross-correlation amplitude, first rank
- X2_AB - cross-correlation amplitude, second rank
- Note: this function is not sensitive to the isotropic components
- of A and B tensors.
- Check consistency
- Components
