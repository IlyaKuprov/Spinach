# kernel/derivatives/fdweights.m

Direct source: [kernel/derivatives/fdweights.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/fdweights.m)
Spin Dynamics Wiki: [fdweights.m](https://spindynamics.org/wiki/index.php?title=fdweights.m)

## Purpose and interface

w=fdweights(target_point,grid_points,max_order) computes finite-difference coefficient rows at a target coordinate from supplied grid coordinates. It uses Fornberg's recursive coefficient construction. Derivative order zero is included, so the first row gives interpolation weights.

- target_point is one real numeric scalar inside the interval from the smallest to largest grid coordinate, including the end points.
- grid_points is a real numeric vector sorted in ascending order.
- max_order is the highest derivative order requested and is less than the number of grid points.
- w has max_order+1 rows and one column per grid point. Row 1 corresponds to order 0; row r+1 supplies order r.

## Recurrence and units

The recurrence begins with the order-zero coefficient for the first grid point and adds grid points one at a time. For each new point it updates the coefficients for the derivative orders available at that stage, through max_order. The resulting row, dotted with function values at grid_points, approximates the requested derivative at target_point.

The grid coordinates determine the scale: for coordinates expressed in units of length, the row for derivative order r has units of inverse length to the r power when applied to function values. No separate spacing parameter or post-scaling is introduced.

## Guards and example

The source requires real numeric input arguments, a scalar target, a vector of grid points, a target lying within their minimum-to-maximum interval, sorted ascending grid points, integer max_order, and max_order < numel(grid_points). The implementation's explicit sortedness test does not separately require strictly distinct grid coordinates. The guard block does not contain separate scalar or nonnegative tests for max_order; it checks integrality and the upper bound only.

    w=fdweights(0,[-1 0 1],1);

Here w has two rows: interpolation weights [0 1 0] and first-derivative weights [-1/2 0 1/2] for these unit-spaced coordinates.

## Related routines

- [fdmat.m](fdmat.md) uses these coefficients to assemble a sparse differentiation matrix.
- [fdvec.m](fdvec.md) uses them to differentiate an input vector.
