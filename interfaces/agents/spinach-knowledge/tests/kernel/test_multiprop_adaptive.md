# tests/kernel/test_multiprop_adaptive.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_multiprop_adaptive.m`
- Signature: `result=test_multiprop_adaptive()`
- Total lines: 140

## Purpose

Tests adaptive repeated propagator application. Syntax: result=test_multiprop_adaptive()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `multiprop()`, `uint64()`, `double()`, `test_close()`, `spdiags()`, `pauli()`, `speye()`, `test_true()`.
