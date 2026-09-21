# tests/kernel/test_dynamic_fp_contexts.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_fp_contexts.m`
- Signature: `result=test_dynamic_fp_contexts()`
- Total lines: 189

## Purpose

Tests compact imaging() and meshflow() context hand-off paths. Syntax: result=test_dynamic_fp_contexts()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_test_imaging()`, `local_test_meshflow()`, `local_context_probe()`, `local_test_square()`, `local_two_cell_mesh()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `imaging()`, `meshflow()`, `local_test_imaging()`, `local_test_meshflow()`, `test_spin_system()`, `state()`, `test_close()`, `local_test_square()`, `test_true()`, `local_two_cell_mesh()`, `speye()`, `local_context_probe()`, `all()`.
