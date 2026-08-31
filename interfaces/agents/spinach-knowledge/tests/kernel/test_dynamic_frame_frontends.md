# tests/kernel/test_dynamic_frame_frontends.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_frame_frontends.m`
- Signature: `result=test_dynamic_frame_frontends()`
- Total lines: 197

## Purpose

Tests frame-transformation and averaging front-end kernels. Syntax: result=test_dynamic_frame_frontends()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file also defines local helper function(s): `local_test_carrier()`, `local_test_frqoffset()`, `local_test_rotframe()`, `local_test_average()`, `local_test_orientation()`, `local_sphten_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Dynamic frame and averaging front ends\n')`.
- Lines 20-23: State the dynamic frame target of the test; implemented by `result=new_test_result('kernel/dynamic_frame_frontends', 'Dynamic frame and averaging front ends', 'frame-transformation helpers must preserve exact algebraic limits.')`.
- Lines 25-26: Check carrier Hamiltonian operator-type algebra; implemented by `result=local_test_carrier(result)`.
- Lines 28-29: Check single-channel and duplicate-channel frequency offsets; implemented by `result=local_test_frqoffset(result)`.
- Lines 31-32: Check the zeroth-order rotating-frame transformation; implemented by `result=local_test_rotframe(result)`.
- Lines 34-35: Check average-Hamiltonian and Krylov-Bogolyubov branches; implemented by `result=local_test_average(result)`.
- Lines 37-38: Check explicit non-zero-angle rotational-basis contraction; implemented by `result=local_test_orientation(result)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_frame_frontends', 'Dynamic frame and averaging front ends', 'frame-transformation helpers must preserve exact algebraic limits.')`.

### Local helper functions

- Line 43: `local_test_carrier()` — `function result=local_test_carrier(result)`. Build a heteronuclear Liouville-space system with non-zero carriers
  - Representative operation: `spin_system=local_sphten_system()`.
  - Representative operation: `basefrqs=spin_system.inter.basefrqs`.
- Line 75: `local_test_frqoffset()` — `function result=local_test_frqoffset(result)`. Build the offset test system and start from a zero Hamiltonian
  - Representative operation: `spin_system=local_sphten_system()`.
  - Representative operation: `H0=sparse(size(spin_system.bas.basis,1),size(spin_system.bas.basis,1))`.
- Line 100: `local_test_rotframe()` — `function result=local_test_rotframe(result)`. Build a one-spin Hilbert-space laboratory-frame system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 124: `local_test_average()` — `function result=local_test_average(result)`. Build a quiet spin system for averaging diagnostics
  - Representative operation: `spin_system=local_sphten_system()`.
  - Representative operation: `Hp=sparse(2,2)`.
- Line 147: `local_test_orientation()` — `function result=local_test_orientation(result)`. Build a sparse rank-one and rank-two rotational basis by hand
  - Representative operation: `Q=cell(2,1)`.
  - Representative operation: `Q{1}=cell(3,3)`.
- Line 182: `local_sphten_system()` — `function spin_system=local_sphten_system()`. Build the common heteronuclear spherical-tensor test system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H','13C'}`.

## Outputs

- result -regression test result with explanatory messages
- The test exercises carrier(), frqoffset(), rotframe(), average(), and
- orientation() on compact systems and analytically controlled limiting
- cases.

## Implementation structure

- Tests frame-transformation and averaging front-end kernels. Syntax:
- result=test_dynamic_frame_frontends()
- result -regression test result with explanatory messages
- The test exercises carrier(), frqoffset(), rotframe(), average(), and
- orientation() on compact systems and analytically controlled limiting
- cases.
- Announce the test target
- State the dynamic frame target of the test
- Check carrier Hamiltonian operator-type algebra
- Check single-channel and duplicate-channel frequency offsets
- Check the zeroth-order rotating-frame transformation
- Check average-Hamiltonian and Krylov-Bogolyubov branches

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_test_carrier()`, `local_test_frqoffset()`, `local_test_rotframe()`, `local_test_average()`, `local_test_orientation()`, `local_sphten_system()`, `carrier()`, `test_close()`, `basefrqs()`, `operator()`, `frqoffset()`, `test_spin_system()`, `assume()`, `rotframe()`, `average()`.
