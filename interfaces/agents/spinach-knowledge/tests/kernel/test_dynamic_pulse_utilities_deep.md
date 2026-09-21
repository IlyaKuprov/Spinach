# tests/kernel/test_dynamic_pulse_utilities_deep.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_pulse_utilities_deep.m`
- Signature: `result=test_dynamic_pulse_utilities_deep()`
- Total lines: 219

## Purpose

Tests dynamic pulse utility paths. Syntax: result=test_dynamic_pulse_utilities_deep()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_liouv_system()`, `local_hilb_system()`, `grad_pulse_ref()`, `grad_sandw_ref()`, `local_delete()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks gradient-pulse dynamics, heterodyne filtering, RLC
- response transforms, Bruker pulse-file writing, finite-RF R-sequence
- compilation, waveform basis variants, and pulse-shape variants.

## Implementation structure

- Tests dynamic pulse utility paths. Syntax:
- result=test_dynamic_pulse_utilities_deep()
- result -regression test result with explanatory messages
- The test checks gradient-pulse dynamics, heterodyne filtering, RLC
- response transforms, Bruker pulse-file writing, finite-RF R-sequence
- compilation, waveform basis variants, and pulse-shape variants.
- Announce the test target
- State the pulse-utility target of the test
- Build a one-proton Liouville-space spin system with a carrier frequency
- Compare grad_pulse with a direct small-matrix exponential reference
- Compare grad_sandw with a direct small-matrix exponential reference
- Heterodyne a clean wall-clock carrier into the rotating frame

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_liouv_system()`, `operator()`, `state()`, `grad_pulse()`, `grad_pulse_ref()`, `test_close()`, `test_true()`, `propagator()`, `grad_sandw()`, `grad_sandw_ref()`, `heterodyne()`, `X_het()`, `Y_het()`, `all()`, `restrans()`.
