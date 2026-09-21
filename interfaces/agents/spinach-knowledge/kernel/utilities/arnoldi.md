# kernel/utilities/arnoldi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/arnoldi.m`
- Signature: `[V,H]=arnoldi(Op,v0,niter)`
- Total lines: 96

## Purpose

Arnoldi procedure for the creation of an orthonormal Krylov basis from repeated action by an operator on a vector. The procedure is numerically unstable and must be used with caution. Syntax: [V,H]=arnoldi(Op,v0,niter)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- Op -function handle taking in a column vector
- and returning another column vector
- v0 -starting vector of the Arnoldi process
- nsteps -number of iterations to take; the Krylov
- subspace will be nsteps+1 dimensional

## Outputs

- V -a matrix containing the orthonormal basis vec-
- tors of the Krylov subspace in columns
- H -extended Hessenberg matrix
- If exact Krylov breakdown occurs, V and H are truncated to the
- completed invariant subspace.

## Implementation structure

- Arnoldi procedure for the creation of an orthonormal Krylov basis
- from repeated action by an operator on a vector. The procedure is
- numerically unstable and must be used with caution. Syntax:
- [V,H]=arnoldi(Op,v0,niter)
- Op -function handle taking in a column vector
- and returning another column vector
- v0 -starting vector of the Arnoldi process
- nsteps -number of iterations to take; the Krylov
- subspace will be nsteps+1 dimensional
- V -a matrix containing the orthonormal basis vec-
- tors of the Krylov subspace in columns
- H -extended Hessenberg matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `iscolumn()`, `isscalar()`.
