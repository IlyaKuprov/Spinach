# kernel/conventions/transforms/qform2sph.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/qform2sph.m`
- Signature: `[r0,r1,r2]=qform2sph(A)`
- Total lines: 58

## Purpose

Returns the spherical harmonic expansion coefficients of the following quadratic form: [x y z]*A*[x y z]'/norm([x y z],2)^2 = sum(r_LM*Y_LM)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check consistency; implemented by `grumble(A)`.
- Lines 31-32: Zeroth rank; implemented by `r0=(2/3)*sqrt(pi)*trace(A)`.
- Lines 34-35: First rank from 1 to -1; implemented by `r1(1)=0; r1(2)=0; r1(3)=0`.
- Lines 37-38: Second rank from 2 to -2; implemented by `r2(1)=+sqrt(2*pi/15)*(A(1,1)-A(2,2)-1i*(A(1,2)+A(2,1)))`.

### Key state/data transformations

- Lines 32: computes `r0` using `r0=(2/3)*sqrt(pi)*trace(A)`.
- Lines 35: computes `r1(1)` using `r1(1)=0; r1(2)=0; r1(3)=0`.
- Lines 38: computes `r2(1)` using `r2(1)=+sqrt(2*pi/15)*(A(1,1)-A(2,2)-1i*(A(1,2)+A(2,1)))`.
- Lines 39: computes `r2(2)` using `r2(2)=-sqrt(2*pi/15)*(A(1,3)+A(3,1)-1i*(A(2,3)+A(3,2)))`.
- Lines 40: computes `r2(3)` using `r2(3)=-(2/3)*sqrt(pi/5)*(-2*A(3,3)+A(2,2)+A(1,1))`.
- Lines 41: computes `r2(4)` using `r2(4)=+sqrt(2*pi/15)*(A(1,3)+A(3,1)+1i*(A(2,3)+A(3,2)))`.
- Lines 42: computes `r2(5)` using `r2(5)=+sqrt(2*pi/15)*(A(1,1)-A(2,2)+1i*(A(1,2)+A(2,1)))`.

### Local helper functions

- Line 47: `grumble()` — `function grumble(A)`. Boring things in physics may be defined as those that make sense.
  - Representative operation: `if (~isnumeric(A))||(size(A,1)~=3)||(size(A,2)~=3)|| (~isreal(A))||(~issymmetric(A))`.
  - Representative operation: `(~isreal(A))||(~issymmetric(A))`.

## Syntax

```matlab
[r0,r1,r2]=qform2sph(A)
```

## Parameters / inputs

- A -a symmetric 3x3 matrix
- Output:
- [r0,r1,r2] -coefficients for zero, first,
- and second rank spherical
- harmonics in the order of
- decreasing m index of Ylm

## Implementation structure

- Returns the spherical harmonic expansion coefficients of the
- following quadratic form:
- [x y z]*A*[x y z]'/norm([x y z],2)^2 = sum(r_LM*Y_LM)
- [r0,r1,r2]=qform2sph(A)
- A -a symmetric 3x3 matrix
- Output:
- [r0,r1,r2] -coefficients for zero, first,
- and second rank spherical
- harmonics in the order of
- decreasing m index of Ylm
- Check consistency
- Zeroth rank

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `issymmetric()`.
