# experiments/pseudocon/interpmat.m

- Signature: `P=interpmat(cube_dims,ranges,xyz)`

## Purpose

Builds a sparse linear operator that interpolates a stretched PCS cube at specified nuclear coordinates using tricubic interpolation. This is an internal part of the PCS inverse-problem solver and is not normally called directly.

## Parameters / inputs

- `cube_dims` — three integer grid dimensions ordered [X Y Z]; each must be at least 4.
- `ranges` — six-element real vector [xmin xmax ymin ymax zmin zmax] in Å, with each minimum less than its maximum.
- `xyz` — N-by-3 real coordinate array in Å. Every point must lie within the supplied ranges; extrapolation is not supported.

## Output

- `P` — sparse matrix with one row per coordinate and `prod(cube_dims)` columns. Multiplying `P` by the flattened PCS cube returns the interpolated values at the corresponding rows of `xyz`.

## Method

The routine creates equally spaced grids on the three axes. For each requested coordinate, it selects four stencil nodes per axis, computes zero-order finite-difference interpolation weights with `fdweights`, and forms their tensor product. These weights are assembled into the rows of the sparse matrix `P`; the per-coordinate construction uses `parfor`.

## References

- [Spin Dynamics Wiki: interpmat.m](https://spindynamics.org/wiki/index.php?title=interpmat.m)
