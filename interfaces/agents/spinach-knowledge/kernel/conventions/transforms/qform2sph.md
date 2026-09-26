# kernel/conventions/transforms/qform2sph.m

- Signature: `[r0,r1,r2]=qform2sph(A)`

## Purpose

Returns the spherical harmonic expansion coefficients of the following quadratic form: [x y z]*A*[x y z]'/norm([x y z],2)^2 = sum(r_LM*Y_LM)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

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
