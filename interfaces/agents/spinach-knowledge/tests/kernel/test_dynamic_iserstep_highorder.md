# tests/kernel/test_dynamic_iserstep_highorder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_iserstep_highorder.m`
- Signature: `result=test_dynamic_iserstep_highorder()`
- Total lines: 89

## Purpose

Tests nonlinear high-order iserstep branches. Syntax: result=test_dynamic_iserstep_highorder()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file also defines local helper function(s): `local_hilb_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks zero-step handling, nonlinear generator execution, and
- agreement of high-order Lie and RKMK branches against a refined DP8
- reference on a compact Hilbert-space problem.

## Implementation structure

- Tests nonlinear high-order iserstep branches. Syntax:
- result=test_dynamic_iserstep_highorder()
- result -regression test result with explanatory messages
- The test checks zero-step handling, nonlinear generator execution, and
- agreement of high-order Lie and RKMK branches against a refined DP8
- reference on a compact Hilbert-space problem.
- Announce the test target
- State the nonlinear Lie-step target of the test
- Build a one-proton Hilbert-space spin system
- Build a Hermitian density matrix with non-zero coherences
- Define a mildly nonlinear, non-commuting Hamiltonian field
- Check the explicit zero-time shortcut in LG4A

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_hilb_system()`, `pauli()`, `iserstep()`, `test_close()`, `tolerances()`, `test_spin_system()`.
