# kernel/derivatives/fdvec.m

- Signature: `dx=fdvec(x,npoints,order)`

## Purpose

Computes the specified derivative of a row or column vector using finite differences. Interior elements use a centered stencil; elements near either end use sided stencils with the same number of points.

## Parameters / inputs

- `x` — numeric row or column vector to differentiate.
- `npoints` — number of points in each finite-difference stencil; must be a positive odd integer.
- `order` — derivative order; must be a positive integer smaller than `npoints`.

## Outputs

- `dx` — derivative values in a vector with the same shape as `x`.

## Numerical / algorithmic content

- Coefficients are computed by `fdweights` for unit-spaced sample positions. No spacing parameter is supplied, so the result is a derivative with respect to the sample index.
- The first `(npoints-1)/2` elements use weights evaluated at their positions within the first `npoints` samples. The corresponding elements at the right end use reversed weights with a factor of `(-1)^order`.
- Interior elements use a centered, symmetric `npoints`-point stencil.

## Validation

- `x` must be a numeric vector. The implementation also checks that it has at least three elements, although its error message says “more than three elements.”
- `npoints` and `order` must satisfy the integer, positivity, odd-stencil and derivative-order constraints above. The implementation does not separately check that `npoints` is no larger than the length of `x`.

Source reference: <https://spindynamics.org/wiki/index.php?title=fdvec.m>
