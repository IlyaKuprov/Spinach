# kernel/utilities/svd_shrink.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/svd_shrink.m`
- Signature: `[vec,cov]=svd_shrink(spin_system,rho,tol)`
- Total lines: 64

## Purpose

Generates sets of vector-covector pairs for the parallel implementation of the time propagation algorithm described in [vec,cov]=svd_shrink(spin_system,rho,tol)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- rho -density matrix
- tol -singluar value drop tolerance

## Outputs

- vec -vectors as columns of a matrix
- cov -covectors as columns of a matrix

## Implementation structure

- Generates sets of vector-covector pairs for the parallel
- implementation of the time propagation algorithm described in
- [vec,cov]=svd_shrink(spin_system,rho,tol)
- rho - density matrix
- tol - singluar value drop tolerance
- vec - vectors as columns of a matrix
- cov - covectors as columns of a matrix
- Check consistency
- Run the singular value decomposition
- Get the drop mask
- Update the user
- Eliminate small singular values

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `report()`, `num2str()`, `nnz()`, `vec()`, `cov()`, `isscalar()`.
