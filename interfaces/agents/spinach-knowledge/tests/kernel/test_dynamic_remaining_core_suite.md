# tests/kernel/test_dynamic_remaining_core_suite.m

- Signature: `result=test_dynamic_remaining_core_suite()`

## Purpose

Tests remaining deterministic utility helpers. Syntax: result=test_dynamic_remaining_core_suite()

## Physical / mathematical content

- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.
## Numerical / algorithmic content

- Uses a Schur-complement reference for adiabatic elimination, analytic checks for coupling and line-shape helpers, and triangle quadrature for a Gaussian convolution. It also tests matrix-based isotope replacement, interaction-representation and frozen-column utilities, and small random-rotation diagnostics.
## Outputs

- `result` — regression test result with explanatory messages.
## Implementation structure

- Check adiabatic elimination against the exact Schur-complement term and verify a Clebsch–Gordan coefficient.
- Check report/banner output to a file and payload impounding.
- Test dipolar coupling and proximity matrices, isotope replacement and its coupling effects, the zeroth-order interaction representation, and frozen-column removal.
- Check Lorentzian and Gaussian line shapes, including small-centre and coincident-vertex cases, and Gaussian triangle quadrature.
- Test magnetic-pumping source-column insertion, `sec2kite` pruning, the Sorensen eigenvalue bound, zero-dynamics trajectory stitching, and the initial orientation from a random-rotation walk.