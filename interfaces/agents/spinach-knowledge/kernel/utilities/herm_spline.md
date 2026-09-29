# kernel/utilities/herm_spline.m

## Purpose

Evaluates a cubic Hermite spline on the [0,1] interval, defined by function values and derivatives at the two interval edges, at one or more query points.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/herm_spline.m>

## Behaviour

- Syntax: `y=herm_spline(f0,df0,f1,df1,x)`.
- The spline coefficients (ordered from x^3 down to x^0) are obtained as a fixed 4x4 matrix multiplying the vector `[df0; f0; df1; f1]`, with the matrix rows `[1 2 1 -2; -2 -3 -1 3; 1 0 0 0; 0 1 0 0]`.
- Each query point is evaluated as `c(1)*x^3+c(2)*x^2+c(3)*x+c(4)` in a loop over all entries of `x`.
- If all four edge values/derivatives are scalars and `x` is not a scalar, the scalars are expanded to arrays of the same size as `x`, so the same spline is evaluated at all query points; otherwise, multiple splines are evaluated at their corresponding query points.
- Input consistency is enforced by an internal `grumble` subfunction:
 all inputs must be numeric, real, and finite (error otherwise), and `f0`, `df0`, `f1`, `df1` must each be either a scalar or the same size as `x` (error otherwise).

## Inputs and outputs

Inputs:

- `f0` - function value(s) at the left edge; real scalar or array.
- `df0` - function derivative(s) at the left edge; real scalar or array.
- `f1` - function value(s) at the right edge; real scalar or array.
- `df1` - function derivative(s) at the right edge; real scalar or array.
- `x` - query point(s) inside the [0,1] interval; real scalar or array.

Output:

- `y` - value of the spline(s) at the query point(s), same size as `x`.

## References

- Spinach Dynamics Wiki: <https://spindynamics.org/wiki/index.php?title=herm_spline.m>
