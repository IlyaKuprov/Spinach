# kernel/utilities/xyz2pd.m

## Purpose

Bins a three-dimensional Cartesian point cloud into a user-specified regular grid and returns raw per-cell counts. The implementation does not divide by the number of points or cell volume. Source: [Spinach repository](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/xyz2pd.m).

## Behaviour

- Syntax: `density=xyz2pd(coords,x_range,y_range,z_range,x_npts,y_npts,z_npts)`.
- Grid cell edges are computed with `linspace` along each axis, producing `npts+1` edges per axis from the given range.
- Each coordinate column is binned with `discretize` to obtain per-axis cell indices.
- Points with any `NaN` bin index (i.e., outside the grid) are discarded via a validity mask.
- Linear cell indices are formed with `sub2ind` over the grid dimensions `[x_npts y_npts z_npts]`.
- Per-cell counts are accumulated with `accumarray` and reshaped into a three-dimensional array of size `x_npts`-by-`y_npts`-by-`z_npts`.
- Input validation (internal `grumble` function) errors out when:
  - `coords` is not a real finite numeric matrix with exactly three columns, or is empty;
  - any range is not a two-element real finite vector with strictly increasing entries;
  - any point count is not a real integer scalar greater than 1.

## Inputs and outputs

Inputs:

- `coords` — N-by-3 array of Cartesian coordinates.
- `x_range` — two-element vector `[xmin xmax]` giving the grid extent along x.
- `y_range` — two-element vector `[ymin ymax]` giving the grid extent along y.
- `z_range` — two-element vector `[zmin zmax]` giving the grid extent along z.
- `x_npts` — number of grid points along x.
- `y_npts` — number of grid points along y.
- `z_npts` — number of grid points along z.

Output:

- `density` — three-dimensional array containing the number of points falling into each grid cell.

## References

- Spinach wiki page: <https://spindynamics.org/wiki/index.php?title=xyz2pd.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/xyz2pd.m>
