# tests/kernel/test_wigner_angular_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_wigner_angular_suite.m`
- Signature: `result=test_wigner_angular_suite()`
- Total lines: 62

## Purpose

Tests angular-momentum coefficient and spherical-function helpers. Syntax: result=test_wigner_angular_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks Clebsch-Gordan coefficients, Wigner symbols, Wigner D
- matrices, and spherical harmonics against elementary exact values.

## Implementation structure

- Tests angular-momentum coefficient and spherical-function helpers. Syntax:
- result=test_wigner_angular_suite()
- result -regression test result with explanatory messages
- The test checks Clebsch-Gordan coefficients, Wigner symbols, Wigner D
- matrices, and spherical harmonics against elementary exact values.
- Announce the test target
- State the angular target of the test
- Coupling two spin-half particles gives triplet and singlet M=0 amplitudes of 1/sqrt(2)
- Wigner 3j values follow from the relation to Clebsch-Gordan coefficients
- Wigner D matrices are unitary representations and reduce to identity for zero rotation
- Spherical harmonics have elementary normalised values

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `clebsch_gordan()`, `wigner_3j()`, `elementary()`, `wigner_6j()`, `wigner()`, `spher_harmon()`.
