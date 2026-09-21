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
