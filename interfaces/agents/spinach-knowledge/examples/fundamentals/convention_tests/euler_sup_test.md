# examples/fundamentals/convention_tests/euler_sup_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/euler_sup_test.m`
- Signature: `euler_sup_test()`
- Total lines: 100

## Purpose

Euler angle superposition tests.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Identity composition test; implemented by `ang=euler_sup([0 0 0],[0 0 0])`.
- Lines 13-14: Random stress test; implemented by `for n=1:2000`.
- Lines 16-17: Draw random Euler angles; implemented by `ang_one=8*pi*(rand(1,3)-0.5)`.
- Lines 20-21: Compose through Euler superposition utility; implemented by `ang_cmp=euler_sup(ang_one,ang_two)`.
- Lines 23-24: Compose through direct matrix multiplication; implemented by `dcm_ref=euler2dcm(ang_two)*euler2dcm(ang_one)`.
- Lines 26-27: Compare the two composite matrices; implemented by `if norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3`.
- Lines 33-34: Singular branch stress test A; implemented by `for n=1:500`.
- Lines 36-37: Draw random near-singular rotations; implemented by `ang_one=[8*pi*(rand-0.5) 1e-12*randn() 8*pi*(rand-0.5)]`.
- Lines 53-54: Singular branch stress test B; implemented by `for n=1:500`.
- Lines 56-57: Draw random near-singular rotations; implemented by `ang_one=[8*pi*(rand-0.5) pi+1e-12*randn() 8*pi*(rand-0.5)]`.
- Lines 73-74: Singular branch stress test C; implemented by `for n=1:500`.
- Lines 76-77: Draw same-phase rotations adding to beta=pi; implemented by `alp=2*pi*(rand-0.5)`.
- Lines 96-97: Successful completion message; implemented by `disp('Euler angle superposition test PASSED.')`.

### Control flow inferred from the code

- Line 9: conditional branch on `norm(euler2dcm(ang)-eye(3),1)>1e-12`.
- Line 14: `for` loop over `n=1:2000`.
- Line 27: conditional branch on `norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3`.
- Line 34: `for` loop over `n=1:500`.
- Line 47: conditional branch on `norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3`.
- Line 54: `for` loop over `n=1:500`.
- Line 67: conditional branch on `norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3`.
- Line 74: `for` loop over `n=1:500`.
- Line 90: conditional branch on `norm(euler2dcm(ang_cmp)-dcm_ref,1)>1e-3`.

### Key state/data transformations

- Lines 8: computes `ang` using `ang=euler_sup([0 0 0],[0 0 0])`.
- Lines 17: computes `ang_one` using `ang_one=8*pi*(rand(1,3)-0.5)`.
- Lines 18: computes `ang_two` using `ang_two=8*pi*(rand(1,3)-0.5)`.
- Lines 21: computes `ang_cmp` using `ang_cmp=euler_sup(ang_one,ang_two)`.
- Lines 24: computes `dcm_ref` using `dcm_ref=euler2dcm(ang_two)*euler2dcm(ang_one)`.
- Lines 77: computes `alp` using `alp=2*pi*(rand-0.5)`.
- Lines 78: computes `gam` using `gam=2*pi*(rand-0.5)`.
- Lines 79: computes `bet` using `bet=pi*rand`.

## Implementation structure

- Euler angle superposition tests.
- Identity composition test
- Random stress test
- Draw random Euler angles
- Compose through Euler superposition utility
- Compose through direct matrix multiplication
- Compare the two composite matrices
- Singular branch stress test A
- Draw random near-singular rotations
- Singular branch stress test B
- Singular branch stress test C
- Draw same-phase rotations adding to beta=pi

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler_sup()`, `euler2dcm()`.
