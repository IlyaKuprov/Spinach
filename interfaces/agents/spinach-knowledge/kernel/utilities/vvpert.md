# kernel/utilities/vvpert.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/vvpert.m`
- Signature: `[Ep,G]=vvpert(E0,H1,order)`
- Total lines: 126

## Purpose

Van Vleck perturbation theory, following Shavitt and Redmon, but excluding the quasi-degenerate split. Syntax: [Ep,G]=vvpert(E0,H1,order)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(E0,H1,order)`.
- Lines 43-44: Reciprocal energy differences; implemented by `Q=1./(E0'-E0); Q(logical(speye(size(Q))))=0`.
- Lines 49-50: Energy differences; implemented by `delta=E0-E0.'`.
- Lines 52-53: Baker-Campbell-Hausdorff coefficients; implemented by `bch=1./factorial(0:order)`.
- Lines 55-56: Allocate recursion storage; implemented by `G=cell(order,1); W=cell(order,1)`.
- Lines 59-60: Compute each perturbation order; implemented by `for n=1:order`.
- Lines 62-63: Build higher nested commutators independent of the current generator; implemented by `for m=2:n`.
- Lines 71-72: Build the transformed Hamiltonian term without the H0 commutator; implemented by `if n==1`.
- Lines 81-82: Split the term into the effective Hamiltonian and generator; implemented by `W{n}=diag(diag(K))`.
- Lines 85-86: Store the first nested commutator for future orders; implemented by `C{2,n+1}=delta.*G{n}`.
- Lines 93-94: Sum the energy corrections; implemented by `Ep=E0`.
- Lines 99-100: Sum the generator corrections; implemented by `G=reshape(G,[1 1 numel(G)])`.

### Control flow inferred from the code

- Line 45: conditional branch on `any(~isfinite(Q),'all')`.
- Line 60: `for` loop over `n=1:order`.
- Line 63: `for` loop over `m=2:n`.
- Line 65: `for` loop over `k=(m-1):(n-1)`.
- Line 72: conditional branch on `n==1`.
- Line 77: `for` loop over `m=2:n`.
- Line 87: conditional branch on `n>1`.
- Line 95: `for` loop over `n=1:order`.

### Key state/data transformations

- Lines 44: computes `Q` using `Q=1./(E0'-E0); Q(logical(speye(size(Q))))=0`.
- Lines 50: computes `delta` using `delta=E0-E0.'`.
- Lines 53: computes `bch` using `bch=1./factorial(0:order)`.
- Lines 56: computes `G` using `G=cell(order,1); W=cell(order,1)`.
- Lines 57: computes `C` using `C=cell(order+1,order+1)`.
- Lines 64: computes `K` using `K=zeros(size(H1),'like',H1)`.
- Lines 68: computes `C{m+1,n+1}` using `C{m+1,n+1}=K`.
- Lines 82: computes `W{n}` using `W{n}=diag(diag(K))`.
- Lines 83: computes `G{n}` using `G{n}=Q.*(K-W{n})`.
- Lines 86: computes `C{2,n+1}` using `C{2,n+1}=delta.*G{n}`.
- Lines 94: computes `Ep` using `Ep=E0`.

### Local helper functions

- Line 106: `grumble()` — `function grumble(E0,H1,order)`.
  - Representative operation: `if (~isnumeric(E0))||(~isreal(E0))||(~iscolumn(E0))`.
  - Representative operation: `error('E0 must be a real column vector.')`.

## Parameters / inputs

- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used; numerical
- artefacts appear beyond about 10-12 for typical
- problems

## Outputs

- Ep -eigenvalues of H0+H1 to the specified order,
- a column vector of reals, not necessarily
- sorted in the same way as the input
- G -Van Vleck transformation generator, such that
- expm(G) is a square unitary matrix with eigen-
- vectors in columns, in the same order as the
- eigenvalues in Ep
- Notes: there must be no degeneracies in H0; H1 must be Hermitian,
- the theory only converges when norm(H1,2) is much smaller
- than the smallest energy gap in H0; complexity is cubic
- both in the order and in the matrix dimension.

## Implementation structure

- Van Vleck perturbation theory, following Shavitt and Redmon, but
- excluding the quasi-degenerate split. Syntax:
- [Ep,G]=vvpert(E0,H1,order)
- E0 -eigenvalues of H0, a column vector of real
- numbers
- H1 -perturbation, written in the basis that di-
- agonalises H0
- order -order of perturbation theory to be used; numerical
- artefacts appear beyond about 10-12 for typical
- problems
- Ep -eigenvalues of H0+H1 to the specified order,
- a column vector of reals, not necessarily

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `logical()`, `speye()`, `any()`, `factorial()`, `comm()`, `bch()`, `cell2mat()`, `iscolumn()`, `ishermitian()`, `isscalar()`.
