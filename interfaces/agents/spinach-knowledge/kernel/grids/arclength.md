# kernel/grids/arclength.m

- Signature: `sig=arclength(r1,r2)`

## Purpose

Returns the central angular distance between two directions on the unit sphere. For a unit sphere, this angle is also the great-circle arc length in sphere-radius units.

## Rule and units

After normalising each input, the routine evaluates `atan2(norm(cross(r1,r2)),dot(r1,r2))`. This gives the shorter central angle from zero to pi radians using the cross-product magnitude and dot product. The scalar result `sig` is in radians; the function does not calculate a weighted grid or a physical length for a sphere with a different radius.

## Parameters / inputs

- `r1`, `r2` — real numeric three-element vectors giving Cartesian endpoint directions. Each supplied vector must have a Euclidean norm within `sqrt(eps)` of one; after this check, the function normalises the vectors and uses column form internally.

## Output

- `sig` — scalar angular distance in radians.

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/arclength.m)
<https://spindynamics.org/wiki/index.php?title=arclength.m>
