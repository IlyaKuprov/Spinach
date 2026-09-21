# tests/kernel/test_pulses_propagation_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_pulses_propagation_suite.m`
- Signature: `result=test_pulses_propagation_suite()`
- Total lines: 128

## Purpose

Tests pulse-coordinate and propagation helpers. Syntax: result=test_pulses_propagation_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Outputs

- result -regression test result with explanatory messages
- The test checks RF coordinate Hessian round-trips, Iserles generators,
- Lie-step methods on a constant generator, and R-sequence phase/compiler
- invariants.

## Implementation structure

- Tests pulse-coordinate and propagation helpers. Syntax:
- result=test_pulses_propagation_suite()
- result -regression test result with explanatory messages
- The test checks RF coordinate Hessian round-trips, Iserles generators,
- Lie-step methods on a constant generator, and R-sequence phase/compiler
- invariants.
- Announce the test target
- State the propagation-helper target of the test
- Choose non-zero amplitudes away from the polar singularity
- Use f=sum(r.^2)=sum(x.^2+y.^2), whose Cartesian Hessian is exactly 2I
- Convert polar coordinates, gradients, and Hessians to Cartesian and back
- Check Iserles second-order and fourth-order product quadrature formulae

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `polar2cartesian()`, `cartesian2polar()`, `test_close()`, `isergen()`, `test_spin_system()`, `operator()`, `state()`, `step()`, `iserstep()`, `rsequence()`, `pauli()`, `rseq_compiler()`, `int2str()`, `speye()`.
