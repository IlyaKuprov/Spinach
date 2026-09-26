# tests/kernel/test_dynamic_iserstep_highorder.m

- Signature: `result=test_dynamic_iserstep_highorder()`

## Purpose

Tests nonlinear high-order iserstep branches. Syntax: result=test_dynamic_iserstep_highorder()

## Physical / mathematical content

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

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
