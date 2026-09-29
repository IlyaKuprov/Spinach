# experiments/pseudocon/interpmat.m

- Signature: `P=interpmat(cube_dims,ranges,xyz)`
- MATLAB source: [`experiments/pseudocon/interpmat.m`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/interpmat.m)

## Purpose

Builds a sparse linear operator that interpolates a stretched PCS density cube at specified Cartesian coordinates. The code describes it as tricubic interpolation and identifies it as an internal part of the PCS inverse-problem solver, not normally called directly.

## Inputs and output

- `cube_dims` — three integer grid sizes in [X Y Z] order; each must be at least 4.
- `ranges` — real six-element vector [xmin xmax ymin ymax zmin zmax] in ångströms (Å), with each minimum less than its corresponding maximum.
- `xyz` — N-by-3 array of real Cartesian coordinates in Å where PCS values are requested; every point must lie inside the supplied ranges.
- `P` — sparse N-by-product(cube_dims) interpolation matrix. Acting on the vectorised cube, it returns one interpolated PCS value per row of `xyz`.

## Interpolation construction

The function creates uniformly spaced x, y, and z grids with `linspace`, using the supplied extents and dimensions. For each query point, it chooses a four-point stencil along each axis, shifting the stencil at grid boundaries so all four indices remain within the cube. It evaluates the one-dimensional interpolation weights with `fdweights` at the requested physical coordinate on each subgrid. The tensor product of the three one-dimensional weight vectors supplies that point's coefficients across the cube; the coefficients are assembled into a sparse matrix row. The point-wise work is in a `parfor` loop.

Queries outside any extent are rejected; extrapolation is not supported. The input checks require three integer dimensions of at least four, real six-element ranges with increasing bounds, and real N-by-3 query coordinates.

## References

- [Spin Dynamics Wiki: interpmat.m](https://spindynamics.org/wiki/index.php?title=interpmat.m)
