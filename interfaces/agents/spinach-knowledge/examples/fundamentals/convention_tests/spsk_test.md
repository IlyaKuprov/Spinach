# examples/fundamentals/convention_tests/spsk_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/spsk_test.m`
- Signature: `spsk_test()`
- Total lines: 34

## Purpose

Test of the span-skew interaction convention.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Eigenvalues; implemented by `xx=rand()`.
- Lines 12-13: Euler angles; implemented by `alp=2*pi*rand()`.
- Lines 17-18: Manual construction; implemented by `R=euler2dcm(alp,bet,gam)`.
- Lines 21-22: Spinach construction; implemented by `iso=(xx+yy+zz)/3`.
- Lines 26-27: Difference; implemented by `if norm(AM-AS,1)<1e-6`.

### Control flow inferred from the code

- Line 27: conditional branch on `norm(AM-AS,1)<1e-6`.

### Key state/data transformations

- Lines 8: computes `xx` using `xx=rand()`.
- Lines 9: computes `yy` using `yy=rand()+3`.
- Lines 10: computes `zz` using `zz=rand()+6`.
- Lines 13: computes `alp` using `alp=2*pi*rand()`.
- Lines 14: computes `bet` using `bet=pi*rand()`.
- Lines 15: computes `gam` using `gam=2*pi*rand()`.
- Lines 18: computes `R` using `R=euler2dcm(alp,bet,gam)`.
- Lines 19: computes `AM` using `AM=R*diag([xx yy zz])*R'`.
- Lines 22: computes `iso` using `iso=(xx+yy+zz)/3`.
- Lines 23: computes `sp` using `sp=zz-xx; sk=3*(yy-iso)/sp`.
- Lines 24: computes `AS` using `AS=spsk2mat(iso,sp,sk,alp,bet,gam)`.

## Implementation structure

- Test of the span-skew interaction convention.
- Eigenvalues
- Euler angles
- Manual construction
- Spinach construction
- Difference

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `spsk2mat()`.
