# tests/kernel/test_dynamic_filter_frontends.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_filter_frontends.m`
- Signature: `result=test_dynamic_filter_frontends()`
- Total lines: 196

## Purpose

Tests dynamic state-filter front-end kernels on compact spin systems. Syntax: result=test_dynamic_filter_frontends()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file also defines local helper function(s): `local_test_coherence_correlation()`, `local_test_decoupling()`, `local_test_homospoil()`, `local_test_spinlock()`, `local_heteronuclear_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Dynamic state-filter front ends\n')`.
- Lines 20-23: State the dynamic filter target of the test; implemented by `result=new_test_result('kernel/dynamic_filter_frontends', 'Dynamic state-filter front ends', 'state-selection front ends must keep exactly the intended basis components.…`.
- Lines 25-26: Check analytical coherence and correlation filters; implemented by `result=local_test_coherence_correlation(result)`.
- Lines 28-29: Check analytical decoupling of a coupled spin system; implemented by `result=local_test_decoupling(result)`.
- Lines 31-32: Check homospoil zero-quantum and longitudinal filters; implemented by `result=local_test_homospoil(result)`.
- Lines 34-35: Check analytical spin-lock projection along both transverse axes; implemented by `result=local_test_spinlock(result)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_filter_frontends', 'Dynamic state-filter front ends', 'state-selection front ends must keep exactly the intended basis components.…`.

### Local helper functions

- Line 40: `local_test_coherence_correlation()` — `function result=local_test_coherence_correlation(result)`. Build a heteronuclear spherical-tensor Liouville-space system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H','13C'}`.
- Line 83: `local_test_decoupling()` — `function result=local_test_decoupling(result)`. Build a coupled heteronuclear system and request the NMR Hamiltonian
  - Representative operation: `spin_system=local_heteronuclear_system()`.
  - Representative operation: `spin_system=assume(spin_system,'nmr')`.
- Line 116: `local_test_homospoil()` — `function result=local_test_homospoil(result)`. Build a homonuclear two-spin system with genuine zero-quantum states
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H','1H'}`.
- Line 149: `local_test_spinlock()` — `function result=local_test_spinlock(result)`. Build a one-spin system for analytical transverse locking
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 181: `local_heteronuclear_system()` — `function spin_system=local_heteronuclear_system()`. Build the common heteronuclear spherical-tensor test system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H','13C'}`.

## Outputs

- result -regression test result with explanatory messages
- The test exercises coherence(), correlation(), decouple(), homospoil(),
- and spinlock() on small spherical-tensor Liouville-space systems with
- analytically known surviving state components.

## Implementation structure

- Tests dynamic state-filter front-end kernels on compact spin systems. Syntax:
- result=test_dynamic_filter_frontends()
- result -regression test result with explanatory messages
- The test exercises coherence(), correlation(), decouple(), homospoil(),
- and spinlock() on small spherical-tensor Liouville-space systems with
- analytically known surviving state components.
- Announce the test target
- State the dynamic filter target of the test
- Check analytical coherence and correlation filters
- Check analytical decoupling of a coupled spin system
- Check homospoil zero-quantum and longitudinal filters
- Check analytical spin-lock projection along both transverse axes

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_test_coherence_correlation()`, `local_test_decoupling()`, `local_test_homospoil()`, `local_test_spinlock()`, `test_spin_system()`, `state()`, `coherence()`, `test_close()`, `correlation()`, `local_heteronuclear_system()`, `assume()`, `hamiltonian()`, `decouple()`, `H_obs()`, `homospoil()`.
