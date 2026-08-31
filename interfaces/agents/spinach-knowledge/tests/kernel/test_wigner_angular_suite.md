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


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Angular-momentum coefficient functions\n')`.
- Lines 19-22: State the angular target of the test; implemented by `result=new_test_result('kernel/wigner_angular_suite', 'Angular-momentum coefficient functions', 'angular-momentum helpers must reproduce elementary exact quantum-mechani…`.
- Lines 24-27: Coupling two spin-half particles gives triplet and singlet M=0 amplitudes of 1/sqrt(2); implemented by `result=test_close(result,'clebsch_gordan triplet half-half', clebsch_gordan(1,0,1/2,1/2,1/2,-1/2),1/sqrt(2),1e-14,1e-14, 'the |1,0> triplet contains |alpha beta> with co…`.
- Lines 35-37: Wigner 3j values follow from the relation to Clebsch-Gordan coefficients; implemented by `result=test_close(result,'wigner_3j 1 1 0',wigner_3j(1,0,1,0,0,0),-1/sqrt(3),1e-14,1e-14, 'the elementary (1 1 0; 0 0 0) Wigner 3j symbol is -1/sqrt(3)')`.
- Lines 43-44: Wigner D matrices are unitary representations and reduce to identity for zero rotation; implemented by `D1=wigner(1,0,pi/2,0)`.
- Lines 55-56: Spherical harmonics have elementary normalised values; implemented by `th=[0 pi/2 pi]; ph=[0 pi/3 pi/7]`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/wigner_angular_suite', 'Angular-momentum coefficient functions', 'angular-momentum helpers must reproduce elementary exact quantum-mechani…`.
- Lines 44: computes `D1` using `D1=wigner(1,0,pi/2,0)`.
- Lines 45: computes `D1_ref` using `D1_ref=[1/2 -1/sqrt(2) 1/2; 1/sqrt(2) 0 -1/sqrt(2); 1/2 1/sqrt(2) 1/2]`.
- Lines 47: computes `'the rank-one Wigner matrix at beta` using `'the rank-one Wigner matrix at beta=pi/2 has the Brink-Satchler closed form')`.
- Lines 48: computes `D0` using `D0=wigner(2,0,0,0)`.
- Lines 49: computes `D` using `D=wigner(2,0.2,0.4,0.7)`.
- Lines 56: computes `th` using `th=[0 pi/2 pi]; ph=[0 pi/3 pi/7]`.

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
