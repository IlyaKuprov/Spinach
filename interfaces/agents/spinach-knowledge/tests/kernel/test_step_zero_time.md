# tests/kernel/test_step_zero_time.m

- Signature: `result=test_step_zero_time()`

## Purpose

Tests zero-duration propagation. Syntax: result=test_step_zero_time()

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Outputs

- result -regression test result with explanatory messages
- The test checks the identity limit of the propagator: a zero time step
- must leave the density matrix exactly unchanged.

## Implementation structure

- Tests zero-duration propagation. Syntax:
- result=test_step_zero_time()
- result -regression test result with explanatory messages
- The test checks the identity limit of the propagator: a zero time step
- must leave the density matrix exactly unchanged.
- Announce the test target
- State the propagation target of the test
- Build a one-proton Hilbert-space spin system
- Propagate an arbitrary Hermitian density matrix for zero time
- Check the identity limit
