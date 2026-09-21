# kernel/utilities/blinv.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/blinv.m`
- Signature: `[Lsq,Dsq]=blinv(A)`
- Total lines: 59

## Purpose

Blicharski's relaxation theory invariants, as given by Equations 20-21 in http://doi.org/10.1515/zna-1972-1012, with an error and a typo corrected in Equation 21. Syntax: [Lsq,Dsq]=blinv(A) where A is the interaction matrix. This function is not sensitive to the trace of the matrix. Parameters: A -a real 3x3 matrix

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismatrix()`, `any()`.
