# kernel/conventions/transforms/euler_equiv.m

- Signature: `answer=euler_equiv(eulers_a,eulers_b,tol)`

## Purpose

Tests whether two ZYZ active Euler-angle sets represent the same rotation to within an angular tolerance. Euler angles are not unique, so the function compares their rotations, not the angle triples themselves.

## Parameters / inputs

- `eulers_a`, `eulers_b`: real, finite three-element Euler-angle vectors `[alpha beta gamma]` in radians, using the ZYZ active convention.
- `tol`: finite, non-negative scalar tolerance in radians.

## Output

- `answer`: true if the relative rotation angle is less than or equal to `tol`.

## Method

The function obtains direction-cosine matrices with `euler2dcm`, forms the relative rotation `dcm_b*dcm_a'`, and compares its geodesic angle on SO(3) to `tol`.

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=euler_equiv.m)
