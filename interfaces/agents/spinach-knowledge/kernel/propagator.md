# kernel/propagator.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/propagator.m`
- Signature: `P=propagator(spin_system,L,timestep)`
- Total lines: 244

## Purpose

Calculates exponential propagator exp(-1i*L*timestep) using scaled and squared Taylor series method. Syntax: P=propagator(spin_system,L,timestep)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(L,timestep); tic`.
- Lines 39-40: Fast bypass for small matrices; implemented by `if size(L,1)<spin_system.tols.small_matrix`.
- Lines 44-45: Load the cache record if one exists; implemented by `if ismember('prop_cache',spin_system.sys.enable)`.
- Lines 47-48: Hash the generator, the step, and the tolerance; implemented by `prop_hash=md5_hash({L,timestep,spin_system.tols.prop_chop})`.
- Lines 51-52: Get ValueStore; implemented by `if ~isworkernode`.
- Lines 58-59: Try to retrieve the propagator from the ValueStore; implemented by `if isKey(store,prop_hash), P=store(prop_hash); return; end`.
- Lines 63-64: Set and clean up a shorthand for -i*L*dt; implemented by `A=clean_up(spin_system,(-1i*timestep)*L,spin_system.tols.prop_chop)`.
- Lines 66-70: Inform the user about generator densities; implemented by `report(spin_system,['generator dimension ' num2str(size(A,1)) ', nnz ' num2str(nnz(A)) ', density ' num2str(100*nnz(A)/numel(A)) '%, sparsity ' num2str(issparse(A))])`.
- Lines 72-73: Estimate the norm; implemented by `mat_norm=cheap_norm(A)`.
- Lines 75-76: Check the norm; implemented by `if mat_norm>1e9`.
- Lines 78-79: The user is doing something silly, bomb out; implemented by `error('norm of -i*L*dt exceeds 1e9, check your L and your dt.')`.
- Lines 83-84: The user is really pushing it, take precautionary measures; implemented by `report(spin_system,'WARNING - time step greatly exceeds the timescale of system dynamics.')`.
- Lines 90-91: Inform the user just in case; implemented by `report(spin_system,'WARNING - time step exceeds the timescale of system dynamics.')`.
- Lines 95-96: Determine scaling and squaring parameters; implemented by `n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings`.
- Lines 100-101: Scale the matrix; implemented by `if scaling_factor>1, A=(1/scaling_factor)*A; end`.
- Lines 103-104: Get the propagator; implemented by `if ismember('gpu',spin_system.sys.enable)&&(size(A,1)>500)`.
- Lines 106-107: Run Taylor series procedure on the GPU; implemented by `A=gpuArray(A); P=speye(size(A))`.
- Lines 111-112: Compute the next term; implemented by `if issparse(A)`.

### Control flow inferred from the code

- Line 40: conditional branch on `size(L,1)<spin_system.tols.small_matrix`.
- Line 45: conditional branch on `ismember('prop_cache',spin_system.sys.enable)`.
- Line 52: conditional branch on `~isworkernode`.
- Line 59: conditional branch on `isKey(store,prop_hash), P=store(prop_hash); return; end`.
- Line 76: conditional branch on `mat_norm>1e9`.
- Line 101: conditional branch on `scaling_factor>1, A=(1/scaling_factor)*A; end`.
- Line 104: conditional branch on `ismember('gpu',spin_system.sys.enable)&&(size(A,1)>500)`.
- Line 109: `while` loop over `nnz(next_term)>0`.
- Line 112: conditional branch on `issparse(A)`.
- Line 133: `while` loop over `nnz(next_term)>0`.
- Line 136: conditional branch on `issparse(A)`.
- Line 168: conditional branch on `n_squarings>0`.
- Line 171: conditional branch on `ismember('gpu',spin_system.sys.enable)`.
- Line 177: `for` loop over `n=1:n_squarings`.

### Key state/data transformations

- Lines 41: computes `P` using `P=expm((-1i*timestep)*L); return`.
- Lines 48: computes `prop_hash` using `prop_hash=md5_hash({L,timestep,spin_system.tols.prop_chop})`.
- Lines 53: computes `store` using `store=gcp('nocreate').ValueStore`.
- Lines 64: computes `A` using `A=clean_up(spin_system,(-1i*timestep)*L,spin_system.tols.prop_chop)`.
- Lines 73: computes `mat_norm` using `mat_norm=cheap_norm(A)`.
- Lines 86: computes `spin_system.tols.prop_chop` using `spin_system.tols.prop_chop=eps('double')`.
- Lines 96: computes `n_squarings` using `n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings`.
- Lines 108: computes `next_term` using `next_term=gpuArray.speye(size(A)); n=1`.

### Local helper functions

- Line 226: `grumble()` — `function grumble(L,timestep)`.
  - Representative operation: `if isa(L,'polyadic')`.
  - Representative operation: `error('L is a polyadic - use evolution() instead.')`.

## Parameters / inputs

- L -Hamiltonian or Liouvillian matrix. If L is
- assembled manually from Hamiltonian commutation
- superoperator H, relaxation superoperator R,
- and kinetics superoperator K, use L=H+1i*R+1i*K
- timestep -propagation time step

## Outputs

- P -propagator matrix
- Note: GPUs are supported, add 'gpu' to sys.enable array during
- calculation setup.
- Note: propagator caching (https://doi.org/10.1063/1.4928978) is
- supported, add 'prop_cache' to sys.enable array to enable.
- Note: we did have Chebyshev and Newton series here at one point,
- as well as the Pade method. None of them had lived up to
- their marketing.

## Implementation structure

- Calculates exponential propagator exp(-1i*L*timestep) using scaled
- and squared Taylor series method. Syntax:
- P=propagator(spin_system,L,timestep)
- L - Hamiltonian or Liouvillian matrix. If L is
- assembled manually from Hamiltonian commutation
- superoperator H, relaxation superoperator R,
- and kinetics superoperator K, use L=H+1i*R+1i*K
- timestep - propagation time step
- P - propagator matrix
- Note: GPUs are supported, add 'gpu' to sys.enable array during
- calculation setup.
- Note: propagator caching (https://doi.org/10.1063/1.4928978) is

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ismember()`, `md5_hash()`, `gcp()`, `getCurrentValueStore()`, `isKey()`, `store()`, `clean_up()`, `report()`, `num2str()`, `nnz()`, `issparse()`, `cheap_norm()`, `eps()`, `log2()`, `gpuArray()`.
