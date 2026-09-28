# kernel/utilities/xyz2pd.m

- Signature: `density=xyz2pd(coords,x_range,y_range,z_range,...`

## Purpose

Bins a three-dimensional Cartesian point cloud onto a user-specified regular grid. The full call is `density=xyz2pd(coords,x_range,y_range,z_range,x_npts,y_npts,z_npts)`.

## Parameters / inputs

- `coords` — an N-by-3 array of real, finite Cartesian coordinates containing at least one point.
- `x_range` — `[xmin xmax]`, the grid extent along the x axis, with `xmin < xmax`.
- `y_range` — `[ymin ymax]`, the grid extent along the y axis, with `ymin < ymax`.
- `z_range` — `[zmin zmax]`, the grid extent along the z axis, with `zmin < zmax`.
- `x_npts`, `y_npts`, `z_npts` — the numbers of grid cells along the respective axes; each must be an integer greater than 1.

## Numerical / algorithmic content

The function constructs equally spaced cell edges from each range using the corresponding point count plus one. It assigns coordinates to cells with `discretize` and discards points outside the grid. It then counts points by cell using `sub2ind` and `accumarray` and reshapes the result to `[x_npts y_npts z_npts]`.

## Outputs

- `density` — a three-dimensional array of the number of points falling into each grid cell. Despite the function description's use of “probability density,” the output is a count, not a normalized density.

## Source

- <https://spindynamics.org/wiki/index.php?title=xyz2pd.m>