# tests/kernel/test_dynamic_pulse_utilities_deep.m

- Signature: `result=test_dynamic_pulse_utilities_deep()`

## Purpose

Tests dynamic pulse utility paths. Syntax: result=test_dynamic_pulse_utilities_deep()

## Physical / mathematical content

- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
