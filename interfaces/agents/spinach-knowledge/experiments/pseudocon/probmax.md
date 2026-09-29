# experiments/pseudocon/probmax.m

- Source: [experiments/pseudocon/probmax.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/probmax.m)
- Wiki: [probmax.m](https://spindynamics.org/wiki/index.php?title=probmax.m)
- Signature: `[x,y,z]=probmax(probden,ranges)`

## Purpose

Returns the coordinate of the largest value in a three-dimensional sampled probability-density array. It locates a sampled maximum; it does not fit or interpolate a continuous maximum.

## Inputs and coordinate mapping

- `probden` is a real numeric 3-D array with dimensions ordered `[X Y Z]`.
- `ranges` is a real six-element vector `[xmin xmax ymin ymax zmin zmax]`. Each lower bound must be strictly less than its matching upper bound.

For each array dimension, the routine builds an inclusive coordinate vector with `linspace` from the corresponding lower and upper bounds, using that dimension's array length. `ndgrid` combines these vectors into coordinate arrays aligned with `probden`.

## Output and constraints

The routine applies MATLAB's linear `max` to `probden(:)` and returns the `x`, `y`, and `z` coordinates at that linear index. If several entries share the maximum, MATLAB's first-maximum behaviour selects the first in linear array order. The validator checks that `ranges` is real numeric with six elements and strictly increasing bounds, and that `probden` is real numeric and three-dimensional; it does not check that values are nonnegative or normalised as a probability density.
