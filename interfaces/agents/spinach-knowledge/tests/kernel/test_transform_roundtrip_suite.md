# tests/kernel/test_transform_roundtrip_suite.m

- Signature: `result=test_transform_roundtrip_suite()`

## Purpose

Tests deterministic coordinate and tensor transforms. Syntax: result=test_transform_roundtrip_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks transformation functions by using exact geometrical
- identities, algebraic inverses, and known tensor decompositions.

## Implementation structure

- Tests deterministic coordinate and tensor transforms. Syntax:
- result=test_transform_roundtrip_suite()
- result -regression test result with explanatory messages
- The test checks transformation functions by using exact geometrical
- identities, algebraic inverses, and known tensor decompositions.
- Announce the test target
- State the transform target of the test
- Direction-cosine matrices must be orthogonal proper rotations
- Quaternion and angle-axis representations must describe the same rotation
- Euler conversion is ill-conditioned in angles, but DCM reconstruction is unique
- Axiality/rhombicity to matrix with zero Euler angles gives the Mehring-order eigenvalues
- Cartesian and irreducible spherical tensor representations are algebraic inverses
