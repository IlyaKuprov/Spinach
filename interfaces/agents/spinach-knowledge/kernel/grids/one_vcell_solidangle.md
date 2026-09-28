# kernel/grids/one_vcell_solidangle.m

- Signature: `S=one_vcell_solidangle(v,centre)`

## Purpose

Computes the solid angle, in radians, of a convex spherical polygon using the method described at https://doi.org/10.1109/TBME.1983.325207.

## Parameters / inputs

- `v`: `3 x n` matrix whose columns are unit vectors giving the polygon vertices.
- `centre`: optional unit-vector coordinates of a centre vertex, supplied as a three-element column vector.

## Output

- `S`: solid angle in radians.

## Numerical / algorithmic content

The function triangulates the polygon using either the first vertex or, when supplied, `centre`. For each triangle `T`, it calculates `atan2(det(T), 1 + sum of the cyclic pairwise column dot products)` and returns twice the sum of those values. When `centre` is supplied, the vertex sequence is closed by appending its first column.

Inputs must be finite, real, correctly shaped, and unit length within `1e-6`.

## Source link

https://spindynamics.org/wiki/index.php?title=one_vcell_solidangle.m