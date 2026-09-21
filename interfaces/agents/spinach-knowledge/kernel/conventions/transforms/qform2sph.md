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
