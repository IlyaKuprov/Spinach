# tests/kernel/test_dynamic_remaining_core_suite.m

- Signature: `result=test_dynamic_remaining_core_suite()`

## Purpose

Tests remaining deterministic utility helpers. Syntax: result=test_dynamic_remaining_core_suite()

## Physical / mathematical content

- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
