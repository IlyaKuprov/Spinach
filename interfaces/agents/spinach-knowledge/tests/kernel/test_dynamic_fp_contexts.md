# tests/kernel/test_dynamic_fp_contexts.m

- Signature: `result=test_dynamic_fp_contexts()`

## Purpose

Tests compact imaging() and meshflow() context hand-off paths. Syntax: result=test_dynamic_fp_contexts()

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Outputs

- result -regression test result with explanatory messages
- The test checks that both contexts assemble finite, correctly sized
- generators and phantom-derived initial and detection states.

## Implementation structure

- Tests compact imaging() and meshflow() context hand-off paths. Syntax:
- result=test_dynamic_fp_contexts()
- result -regression test result with explanatory messages
- The test checks that both contexts assemble finite, correctly sized
- generators and phantom-derived initial and detection states.
- Announce the test target
- State the context target of the test
- Exercise the Cartesian-grid imaging context
- Exercise the unstructured-mesh flow context
- Build a one-spin spherical-tensor Liouville-space system
- Set a minimal one-dimensional imaging grid
- Supply relaxation, initial-state, and coil phantoms
