# tests/kernel/test_dynamic_remaining_core_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_remaining_core_suite.m`
- Signature: `result=test_dynamic_remaining_core_suite()`
- Total lines: 287

## Purpose

Tests remaining deterministic utility helpers. Syntax: result=test_dynamic_remaining_core_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_liouvillian_system()`, `local_dipolar_system()`, `local_isoswap_inputs()`, `gcp()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks block eliminations, text reporting, spin metadata,
- analytical line shapes, pumping terms, kite pruning, trajectory
- stitching, and small random-rotation diagnostics.

## Implementation structure

- Tests remaining deterministic utility helpers. Syntax:
- result=test_dynamic_remaining_core_suite()
- result -regression test result with explanatory messages
- The test checks block eliminations, text reporting, spin metadata,
- analytical line shapes, pumping terms, kite pruning, trajectory
- stitching, and small random-rotation diagnostics.
- Announce the test target
- State the utility target of the test
- Build a minimal Liouville-space descriptor for algebraic utilities
- Check adiabatic elimination against the exact Schur complement term
- Check a textbook Clebsch-Gordan coefficient
- Check direct console and banner reporting into a file handle

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_liouvillian_system()`, `adelim()`, `test_close()`, `cg_fast()`, `fopen()`, `report()`, `banner()`, `fclose()`, `fileread()`, `delete()`, `test_true()`, `contains()`, `polinfo()`, `polyadic()`, `speye()`.
