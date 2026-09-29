# kernel/grids/sphtrsubd.m

- Signature: `[r12,r23,r31]=sphtrsubd(r1,r2,r3)`

## Purpose

Computes the three great-circle arc midpoints of a spherical triangle.

## Inputs and outputs

- `r1`, `r2`, and `r3` are real, three-element unit vectors giving the triangle vertices.
- `r12`, `r23`, and `r31` are three-element unit vectors for the midpoints of arcs 1–2, 2–3, and 3–1, respectively.

## Construction and guards

For each pair, the function adds the endpoint vectors and normalises the sum to unit length. This is the midpoint on the shorter great-circle arc. Before construction, the helper requires each input to be numeric, real, and have three elements; each norm must differ from 1 by no more than `sqrt(eps)`. It rejects a parent spherical-triangle area above `pi/2` and any of its three arc lengths above `pi/2`. The side-length guard keeps each midpoint sum away from the zero vector.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/sphtrsubd.m)
- [Spinach Wiki: sphtrsubd.m](https://spindynamics.org/wiki/index.php?title=sphtrsubd.m)
