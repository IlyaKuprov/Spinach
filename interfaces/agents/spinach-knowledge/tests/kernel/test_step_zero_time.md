# tests/kernel/test_step_zero_time.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_step_zero_time.m`
- Signature: `result=test_step_zero_time()`
- Total lines: 43

## Purpose

Tests zero-duration propagation. Syntax: result=test_step_zero_time()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `pauli()`, `step()`, `test_close()`, `rho()`.
