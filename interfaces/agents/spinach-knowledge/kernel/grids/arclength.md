# kernel/grids/arclength.m

- Signature: `sig=arclength(r1,r2)`

## Purpose

Returns the arc length between two points on the unit sphere represented by unit vectors.

## Physical / mathematical content

- The function returns the angular distance in radians, computed as `atan2(norm(cross(r1,r2)),dot(r1,r2))`.
- Inputs must be real numeric three-element vectors whose Euclidean norms differ from one by no more than `sqrt(eps)`; the vectors are then normalized before the angle is evaluated.

## Parameters / inputs

- `r1`, `r2` — three-element unit vectors giving the Cartesian coordinates of the arc endpoints.

## Outputs

- `sig` — arc length in radians.

<https://spindynamics.org/wiki/index.php?title=arclength.m>
