# kernel/conventions/transforms/qter2anax.m

- Signature: `[rot_axis,rot_angle]=qter2anax(q)`

## Purpose

Converts a quaternion representing a rotation to its angle and axis.

## Physical / mathematical content

For the normalized quaternion, let `v=[q.i q.j q.k]` and `n=norm(v,2)`. If n is nonzero, `rot_angle=2*atan2(n,q.u)` in radians and `rot_axis=v/n`. If n is zero, the function returns angle 0 and axis [0 0 1].

## Numerical / algorithmic content

The quaternion is normalized by the Euclidean norm of its four components before the angle-axis calculation. Missing fields or non-real quaternion components cause an error.

## Syntax

```matlab
[rot_axis,rot_angle]=qter2anax(q)
```

## Parameters / inputs

- `q` — quaternion structure with real-valued scalar fields `u`, `i`, `j`, and `k`.

## Outputs

- `rot_axis` — three-element row vector giving the rotation axis.
- `rot_angle` — rotation angle in radians.

## Implementation structure

The function normalizes q, computes the norm of its vector part, and handles the zero-vector case separately before returning the angle and axis.
