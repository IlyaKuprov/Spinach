# examples/fundamentals/convention_tests/rotations_1.m

- Signature: `rotations_1()`

## Purpose

Tests the internal consistency of kernel rotation functions.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

## Implementation structure

- Tests the internal consistency of kernel rotation functions.
- Generate a random symmetric traceless 3x3 matrix
- Generate a random set of Euler angles
- % Test 1: euler2dcm, wigner, mat2sphten
- DCM rotation followed by a transformation into irreducible components
- Transformation into irreducible components followed by a Wigner rotation
- Check the difference
- % Test 2: euler2dcm, dcm2euler
- Transforms Euler angles into DCM
- Transform the DCM back into Euler angles
- % Test 3: euler2dcm, wigner, dcm2wigner
- Transforms Euler angles into DCM, then DCM to Wigner matrix
