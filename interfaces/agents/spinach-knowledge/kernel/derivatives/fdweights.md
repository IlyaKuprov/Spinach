# kernel/derivatives/fdweights.m

- Signature: `w=fdweights(target_point,grid_points,max_order)`

## Purpose

Calculates finite-difference weights for numerical derivatives at `target_point` using values at `grid_points`. Derivative order 0 gives interpolation weights.

## Parameters / inputs

- `target_point`: point at which the derivative is required.
- `grid_points`: points at which the function is given, sorted in ascending order.
- `max_order`: maximum derivative order; an integer smaller than the number of grid points.

## Output

- `w`: finite-difference coefficient array with one column per grid point. Row `k+1` contains the weights for derivative order `k`, from 0 through `max_order`.

## Validation and algorithm

The inputs must be real and numeric. `target_point` must be a scalar within the range of `grid_points`, and `grid_points` must be a vector sorted in ascending order. The routine initializes the order-0 weight for the first point, then adds grid points one at a time, updating the weights for all derivative orders up to `max_order`.

[Source reference](https://spindynamics.org/wiki/index.php?title=fdweights.m)