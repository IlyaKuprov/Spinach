# tests/kernel/test_dynamic_propagation_frontends.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_propagation_frontends.m`
- Signature: `result=test_dynamic_propagation_frontends()`
- Total lines: 246

## Purpose

Tests dynamic propagation front-end kernels on tiny systems. Syntax: result=test_dynamic_propagation_frontends()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_test_propagator_step()`, `local_test_evolution()`, `local_test_krylov()`, `local_test_reduce()`, `local_liouville_system()`, `local_hilbert_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Dynamic propagation front ends\n')`.
- Lines 19-22: State the dynamic propagation target of the test; implemented by `result=new_test_result('kernel/dynamic_propagation_frontends', 'Dynamic propagation front ends', 'propagation front ends must match direct finite-dimensional references.…`.
- Lines 24-25: Check scaled Taylor propagator and step() branches; implemented by `result=local_test_propagator_step(result)`.
- Lines 27-28: Check evolution() output modes against explicit propagator products; implemented by `result=local_test_evolution(result)`.
- Lines 30-31: Check direct krylov() output modes against step() references; implemented by `result=local_test_krylov(result)`.
- Lines 33-34: Check reduce() projector invariants and blanket-disable branch; implemented by `result=local_test_reduce(result)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_propagation_frontends', 'Dynamic propagation front ends', 'propagation front ends must match direct finite-dimensional references.…`.

### Local helper functions

- Line 39: `local_test_propagator_step()` — `function result=local_test_propagator_step(result)`. Build a one-spin Liouville-space system and force Taylor propagation
  - Representative operation: `spin_system=local_liouville_system()`.
  - Representative operation: `spin_system.tols.small_matrix=2`.
- Line 93: `local_test_evolution()` — `function result=local_test_evolution(result)`. Build a one-spin Liouville-space system with trajectory reduction disabled
  - Representative operation: `spin_system=local_liouville_system()`.
  - Representative operation: `spin_system.sys.disable=unique([spin_system.sys.disable {'trajlevel'}])`.
- Line 146: `local_test_krylov()` — `function result=local_test_krylov(result)`. Build a one-spin Liouville-space system for direct Krylov propagation
  - Representative operation: `spin_system=local_liouville_system()`.
  - Representative operation: `L=operator(spin_system,'Lz','1H')`.
- Line 187: `local_test_reduce()` — `function result=local_test_reduce(result)`. Build a one-spin Liouville-space system for projector checks
  - Representative operation: `spin_system=local_liouville_system()`.
  - Representative operation: `L=operator(spin_system,'Lz','1H')`.
- Line 220: `local_liouville_system()` — `function spin_system=local_liouville_system()`. Build the common one-spin spherical-tensor Liouville-space system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 234: `local_hilbert_system()` — `function spin_system=local_hilbert_system()`. Build the common one-spin Hilbert-space system
  - Representative operation: `sys.magnet=0.0`.
  - Representative operation: `sys.isotopes={'1H'}`.

## Outputs

- result -regression test result with explanatory messages
- The test exercises propagator(), step(), evolution(), krylov(), and
- reduce() against direct finite-dimensional propagation references.

## Implementation structure

- Tests dynamic propagation front-end kernels on tiny systems. Syntax:
- result=test_dynamic_propagation_frontends()
- result -regression test result with explanatory messages
- The test exercises propagator(), step(), evolution(), krylov(), and
- reduce() against direct finite-dimensional propagation references.
- Announce the test target
- State the dynamic propagation target of the test
- Check scaled Taylor propagator and step() branches
- Check evolution() output modes against explicit propagator products
- Check direct krylov() output modes against step() references
- Check reduce() projector invariants and blanket-disable branch
- Build a one-spin Liouville-space system and force Taylor propagation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_test_propagator_step()`, `local_test_evolution()`, `local_test_krylov()`, `local_test_reduce()`, `local_liouville_system()`, `operator()`, `state()`, `propagator()`, `test_close()`, `step()`, `iserstep()`, `local_hilbert_system()`, `pauli()`, `evolution()`, `traj_ref()`.
