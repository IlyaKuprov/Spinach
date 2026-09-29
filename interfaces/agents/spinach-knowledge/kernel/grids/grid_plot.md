# kernel/grids/grid_plot.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/grid_plot.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=grid_plot.m)

## Purpose

Plot spherical sample points and their Voronoi tessellation. The routine draws into the current graphics axes and returns no value.

## Inputs

- `x`, `y`, and `z` are equal-length N-by-1 column vectors of finite real numeric Cartesian coordinates. They represent points on a sphere; the function does not rescale them or check that each point has unit norm, and the coordinates have no physical units assigned by this routine.
- `vorn` is a cell array of tessera vertex data. If omitted or empty, the tessellation is obtained from `voronoisphere([x';y';z'])`. The source checks that `vorn` is a cell array, but does not check that its cell count matches N.
- `c` selects face colour: omitted or empty means white; a character value is used as the face colour name; otherwise each face uses `c(k)`. Numeric colour data should therefore provide one value per face; the routine does not validate its length.
- `options.dots` controls centre markers. If the options argument or this field is absent, dots default to true; setting it false suppresses them.

## Rendering behaviour

When enabled, the sample centres are black dots with marker size 3. Each tessera is drawn as a fully opaque patch. The plot limits are [-1.1,+1.1] on each Cartesian axis, the axes are square, the camera position is [0,0,10], and tick marks are hidden. The routine computes a tessellation only when `vorn` is omitted or empty; it does not calculate quadrature weights.
