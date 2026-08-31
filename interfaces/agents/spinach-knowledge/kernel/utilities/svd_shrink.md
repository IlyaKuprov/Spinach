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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Check consistency; implemented by `grumble(rho,tol)`.
- Lines 29-30: Run the singular value decomposition; implemented by `[vec,S,cov]=svd(full(rho)); S=diag(S)`.
- Lines 32-33: Get the drop mask; implemented by `drop_mask=(S<tol)`.
- Lines 35-37: Update the user; implemented by `report(spin_system,['dropped ' num2str(nnz(drop_mask)) ' insignificant ' 'vector-covector pairs from the density matrix.'])`.
- Lines 39-40: Eliminate small singular values; implemented by `vec(:,drop_mask)=[]; cov(:,drop_mask)=[]; S(drop_mask)=[]`.
- Lines 42-43: Spread the coefficients; implemented by `vec=vec*diag(sqrt(S))`.

### Key state/data transformations

- Lines 30: computes `[vec,S,cov]` using `[vec,S,cov]=svd(full(rho)); S=diag(S)`.
- Lines 33: computes `drop_mask` using `drop_mask=(S<tol)`.
- Lines 40: computes `vec(:,drop_mask)` using `vec(:,drop_mask)=[]; cov(:,drop_mask)=[]; S(drop_mask)=[]`.
- Lines 43: computes `vec` using `vec=vec*diag(sqrt(S))`.
- Lines 44: computes `cov` using `cov=cov*diag(sqrt(S))`.

### Local helper functions

- Line 49: `grumble()` — `function grumble(rho,tol)`. If everyone likes your research, you can be certain that you have not done anything important. That is the first
  - Representative operation: `if (~isnumeric(rho))||(size(rho,1)~=size(rho,2))`.
  - Representative operation: `error('rho must be a square matrix.')`.

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
