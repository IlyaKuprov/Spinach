# tests/kernel/test_multiprop_adaptive.m

- Signature: `result=test_multiprop_adaptive()`

## Purpose

Tests adaptive repeated propagator application. Syntax: result=test_multiprop_adaptive()

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- result -regression test result with explanatory messages
- The test checks binary adaptive squaring in multiprop() against explicit
- matrix-power references for state vectors and density matrices, and
- verifies clean-up after propagator squaring.

## Implementation structure

- Tests adaptive repeated propagator application. Syntax:
- result=test_multiprop_adaptive()
- result -regression test result with explanatory messages
- The test checks binary adaptive squaring in multiprop() against explicit
- matrix-power references for state vectors and density matrices, and
- verifies clean-up after propagator squaring.
- Announce the test target
- State the propagation target of the test
- Build the minimum Spinach system fields required by clean_up()
- Define a non-normal sparse propagator and a state vector
- Compare state-vector propagation with an explicit matrix power
- Define a diagonal sparse propagator and a state vector
