# examples/fundamentals/convention_tests/rotations_1.m

- Signature: `rotations_1()`

## Purpose

Checks consistency among Spinach's DCM, Euler-angle, Wigner-matrix, and Cartesian-to-spherical-tensor rotation routines.

## Checks

A random symmetric traceless 3×3 matrix and random Euler angles are used for three comparisons:

1. Rotating the matrix by its DCM and then converting it with `mat2sphten` is compared with converting first and applying `wigner(2,...)`; the rank-2 coefficient residual must be below 10⁻¹⁰ in the 2-norm.
2. Converting the Euler angles to a DCM and back with `dcm2euler` must reproduce the angles within a 2-norm tolerance of 10⁻³.
3. `dcm2wigner(euler2dcm(...))` is compared with the direct `wigner(2,...)` result, with a 2-norm tolerance of 10⁻¹⁰.
