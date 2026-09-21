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
