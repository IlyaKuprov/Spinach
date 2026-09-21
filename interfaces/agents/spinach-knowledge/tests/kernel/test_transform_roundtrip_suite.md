# tests/kernel/test_transform_roundtrip_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_transform_roundtrip_suite.m`
- Signature: `result=test_transform_roundtrip_suite()`
- Total lines: 85

## Purpose

Tests deterministic coordinate and tensor transforms. Syntax: result=test_transform_roundtrip_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `anax2dcm()`, `test_close()`, `anax2qter()`, `qter2anax()`, `euler2dcm()`, `dcm2euler()`, `axrh2mat()`, `mat2axrh()`, `eigvals()`, `mat2sphten()`, `sphten2mat()`, `frac2cart()`, `xyz2sph()`.
