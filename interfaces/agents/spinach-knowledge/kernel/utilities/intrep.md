# kernel/utilities/intrep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/intrep.m`
- Signature: `Hr=intrep(spin_system,H0,H,T,order)`
- Total lines: 113

## Purpose

Interaction representation transformation with respect to a specified Hamiltonian to specified order in perturbation theory (https://doi.org/10.1063/1.4928978). Syntax: Hr=intrep(spin_system,H0,H,T,order)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(H0,H,order)`.
- Lines 38-39: Confirm that T is indeed the period; implemented by `P=propagator(spin_system,H0,T)`.
- Lines 44-45: Compute the norms; implemented by `norm_h0=norm(H0,1); norm_h1=norm(H-H0,1)`.
- Lines 47-48: Report the norms; implemented by `report(spin_system,['rotating frame period: ' num2str(T) ' seconds'])`.
- Lines 52-53: Confirm that norms are sensible; implemented by `if norm_h1>norm_h0`.
- Lines 57-58: Decide the theory order; implemented by `switch order`.
- Lines 62-63: Shortcut for high field; implemented by `Hr=H-H0`.
- Lines 67-68: Shortcut for infinite order; implemented by `Hr=(1i/T)*logm(expm(full(-1i*H*T)))`.
- Lines 72-73: Get the derivatives; implemented by `D=dirdiff(spin_system,H0,H-H0,T,order+1)`.
- Lines 75-76: Get the first term; implemented by `Hr=(1i/T)*(D{1}'*D{2})`.
- Lines 78-79: Get the rest of the series; implemented by `for n=2:order`.
- Lines 87-88: Clean up the output; implemented by `Hr=clean_up(spin_system,(Hr+Hr')/2,spin_system.tols.liouv_zero)`.
- Lines 90-93: Return matrix density statistics; implemented by `report(spin_system,['int. rep. Hamiltonian dimension ' num2str(size(Hr,1)) ', nnz ' num2str(nnz(Hr)) ', density ' num2str(100*nnz(Hr)/numel(Hr)) '%, sparsity ' num2str(i…`.

### Control flow inferred from the code

- Line 40: conditional branch on `norm(P-speye(size(P)),1)>1e-6`.
- Line 53: conditional branch on `norm_h1>norm_h0`.
- Line 58: dispatches on `order`; cases `0`, `Inf`.
- Line 79: `for` loop over `n=2:order`.
- Line 80: `for` loop over `k=1:n`.

### Key state/data transformations

- Lines 39: computes `P` using `P=propagator(spin_system,H0,T)`.
- Lines 45: computes `norm_h0` using `norm_h0=norm(H0,1); norm_h1=norm(H-H0,1)`.
- Lines 63: computes `Hr` using `Hr=H-H0`.
- Lines 73: computes `D` using `D=dirdiff(spin_system,H0,H-H0,T,order+1)`.

### Local helper functions

- Line 98: `grumble()` — `function grumble(H0,H,order)`. There is no difference between communism and socialism, except in the means of achieving the same ultimate end: communism proposes to enslave
  - Representative operation: `if (~ishermitian(H))||(~ishermitian(H0))`.
  - Representative operation: `error('both H and H0 must be Hermitian.')`.

## Parameters / inputs

- H0 -the Hamiltonian with respect to which the
- interaction representation transformation
- is to be done, typically Zeeman Hamiltonian
- H -laboratory frame Hamiltonian H0+H1 that is
- to be transformed into the interaction rep-
- resentation, typically the full Hamiltonian
- T -period of the H0 propagator
- order -perturbation theory order in the rotating
- frame transformation, this may be inf

## Outputs

- Hr -Hamiltonian in the interaction representation
- Note: the auxiliary matrix method is massively faster than
- either commutator series or diagonalisation.

## Implementation structure

- Interaction representation transformation with respect to
- a specified Hamiltonian to specified order in perturbation
- theory (https://doi.org/10.1063/1.4928978). Syntax:
- Hr=intrep(spin_system,H0,H,T,order)
- H0 -the Hamiltonian with respect to which the
- interaction representation transformation
- is to be done, typically Zeeman Hamiltonian
- H -laboratory frame Hamiltonian H0+H1 that is
- to be transformed into the interaction rep-
- resentation, typically the full Hamiltonian
- T -period of the H0 propagator
- order -perturbation theory order in the rotating

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `propagator()`, `speye()`, `report()`, `num2str()`, `logm()`, `dirdiff()`, `nchoosek()`, `factorial()`, `clean_up()`, `nnz()`, `issparse()`, `ishermitian()`, `isinf()`.
