# kernel/grids/sphtrsubd.m

- Signature: `[r12,r23,r31]=sphtrsubd(r1,r2,r3)`

## Purpose

Returns the arc midpoints of a spherical triangle specified by its three vertex unit vectors.

## Inputs

- `r1`, `r2`, `r3`: three-element real unit vectors giving Cartesian coordinates of the triangle vertices.

## Outputs

- `r12`, `r23`, `r31`: three-element unit vectors giving Cartesian coordinates of the arc midpoints for vertex pairs 1–2, 2–3, and 3–1.

## Algorithm and constraints

Each midpoint is computed by adding the corresponding vertex vectors and dividing by the Euclidean norm of the sum. Before calculation, the function checks that each input is a three-element real numeric unit vector, with unit length checked to within `sqrt(eps)`. It rejects parent triangles with area greater than `pi/2` or any vertex-pair arc length greater than `pi/2`.

Source reference: https://spindynamics.org/wiki/index.php?title=sphtrsubd.m