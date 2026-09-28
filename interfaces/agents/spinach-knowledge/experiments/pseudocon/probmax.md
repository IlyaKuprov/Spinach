# experiments/pseudocon/probmax.m

- Signature: `[x,y,z]=probmax(probden,ranges)`

## Purpose

Returns the coordinates of the maximum value in a three-dimensional probability-density array, using the coordinate bounds in `ranges`.

## Numerical / algorithmic content

The routine builds one coordinate grid per dimension with `linspace` and `ndgrid`, matching each grid length to the corresponding dimension of `probden`. It locates the maximum array value and returns the coordinates at its linear index.

## Parameters / inputs

- `probden` — real three-dimensional probability-density array, with dimensions ordered `[X Y Z]`.
- `ranges` — six-element vector `[xmin xmax ymin ymax zmin zmax]` defining the coordinate bounds.

## Outputs

- `x`, `y`, `z` — coordinates of the maximum point.
