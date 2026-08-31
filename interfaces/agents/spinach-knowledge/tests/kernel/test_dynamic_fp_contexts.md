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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Fokker-Planck imaging and meshflow contexts\n')`.
- Lines 19-22: State the context target of the test; implemented by `result=new_test_result('kernel/dynamic_fp_contexts', 'Fokker-Planck imaging and meshflow contexts', 'imaging() and meshflow() must assemble context generators and phanto…`.
- Lines 24-25: Exercise the Cartesian-grid imaging context; implemented by `result=local_test_imaging(result)`.
- Lines 27-28: Exercise the unstructured-mesh flow context; implemented by `result=local_test_meshflow(result)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_fp_contexts', 'Fokker-Planck imaging and meshflow contexts', 'imaging() and meshflow() must assemble context generators and phanto…`.

### Local helper functions

- Line 33: `local_test_imaging()` — `function result=local_test_imaging(result)`. Build a one-spin spherical-tensor Liouville-space system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 90: `local_test_meshflow()` — `function result=local_test_meshflow(result)`. Build a one-spin spherical-tensor Liouville-space system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 141: `local_context_probe()` — `function answer=local_context_probe(~,parameters,H,R,K,G,F)`. Return context products for regression checks
  - Representative operation: `answer.spc_dim=parameters.spc_dim`.
  - Representative operation: `answer.spn_dim=parameters.spn_dim`.
- Line 158: `local_test_square()` — `function result=local_test_square(result,label,A,matrix_dim)`. Check that a context object is a square matrix of the expected size
  - Representative operation: `result=test_close(result,label,size(A),[matrix_dim matrix_dim],0,0, 'context generators must act on the full Fokker-Planck state vector')`.
  - Representative operation: `'context generators must act on the full Fokker-Planck state vector')`.
- Line 169: `local_two_cell_mesh()` — `function mesh=local_two_cell_mesh()`. Active mesh vertices are the parents of the two finite-volume cells
  - Representative operation: `mesh.idx.active=[1 2]`.
  - Representative operation: `mesh.idx.triangles=[1 2 3]`.

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
