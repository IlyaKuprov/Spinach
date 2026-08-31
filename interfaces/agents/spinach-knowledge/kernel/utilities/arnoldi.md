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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(Op,v0,niter)`.
- Lines 36-37: Preallocate outputs; implemented by `V=zeros(numel(v0),niter+1,'like',1i)`.
- Lines 40-41: First vector; implemented by `V(:,1)=v0/norm(v0,2)`.
- Lines 43-44: Iteration loop; implemented by `for n=1:niter`.
- Lines 46-47: Get the next vector; implemented by `V(:,n+1)=Op(V(:,n))`.
- Lines 49-50: Gram-Schmidt loop; implemented by `for k=1:n`.
- Lines 52-53: Get the scalar product; implemented by `H(k,n)=V(:,k)'*V(:,n+1)`.
- Lines 55-56: Remove its contribution; implemented by `V(:,n+1)=V(:,n+1)-H(k,n)*V(:,k)`.
- Lines 60-61: Compute the norm; implemented by `H(n+1,n)=norm(V(:,n+1),2)`.
- Lines 63-64: Return if exact Krylov breakdown occurs; implemented by `if H(n+1,n)==0`.
- Lines 68-69: Divide out the norm; implemented by `V(:,n+1)=V(:,n+1)/H(n+1,n)`.

### Control flow inferred from the code

- Line 44: `for` loop over `n=1:niter`.
- Line 50: `for` loop over `k=1:n`.
- Line 64: conditional branch on `H(n+1,n)==0`.

### Key state/data transformations

- Lines 37: computes `V` using `V=zeros(numel(v0),niter+1,'like',1i)`.
- Lines 38: computes `H` using `H=zeros(niter+1,niter,'like',1i)`.
- Lines 41: computes `V(:,1)` using `V(:,1)=v0/norm(v0,2)`.
- Lines 47: computes `V(:,n+1)` using `V(:,n+1)=Op(V(:,n))`.
- Lines 53: computes `H(k,n)` using `H(k,n)=V(:,k)'*V(:,n+1)`.
- Lines 61: computes `H(n+1,n)` using `H(n+1,n)=norm(V(:,n+1),2)`.

### Local helper functions

- Line 76: `grumble()` — `function grumble(Op,v0,niter)`.
  - Representative operation: `if ~isa(Op,'function_handle')`.
  - Representative operation: `error('Op must be a function handle.')`.

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
