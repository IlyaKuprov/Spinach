# kernel/utilities/blinv.m

- Signature: `[Lsq,Dsq]=blinv(A)`

## Purpose

Blicharski's relaxation theory invariants, as given by Equations 20-21 in http://doi.org/10.1515/zna-1972-1012, with an error and a typo corrected in Equation 21. Syntax: [Lsq,Dsq]=blinv(A) where A is the interaction matrix. This function is not sensitive to the trace of the matrix. Parameters: A -a real 3x3 matrix

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Outputs

- Lsq -first rank invariant
- Dsq -second rank invariant

## Implementation structure

- Blicharski's relaxation theory invariants, as given by Equations
- 20-21 in http://doi.org/10.1515/zna-1972-1012, with an error and
- a typo corrected in Equation 21. Syntax:
- [Lsq,Dsq]=blinv(A)
- where A is the interaction matrix. This function is not sensitive
- to the trace of the matrix. Parameters:
- A - a real 3x3 matrix
- Lsq - first rank invariant
- Dsq - second rank invariant
- Check consistency
- First rank invariant
- Second rank invariant
