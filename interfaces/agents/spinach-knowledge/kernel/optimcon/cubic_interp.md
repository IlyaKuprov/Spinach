# kernel/optimcon/cubic_interp.m

- Signature: `[alpha,fx]=cubic_interp(end_a,end_b,alpha_a,alpha_b,...`

## Purpose

Finds the maximum of a cubic interpolant within the interval bounded by `end_a` and `end_b`. The interpolant is determined by function values and directional derivatives at `alpha_a` and `alpha_b`.

## Numerical / algorithmic content

The function constructs the cubic in coordinates normalised relative to `alpha_a` and `alpha_b`. It finds the roots of the cubic's derivative, discards complex roots and roots outside the interpolation interval, and compares the cubic values at the remaining roots and both interval boundaries. It returns the candidate with the largest value, transformed back to alpha coordinates.

## Parameters / inputs

- `end_a` — first interpolation boundary in alpha space.
- `end_b` — second interpolation boundary in alpha space.
- `alpha_a` — first interpolation anchor point.
- `alpha_b` — second interpolation anchor point; must differ from `alpha_a`.
- `f_a` — function value at `alpha_a`.
- `dir_der_a` — directional derivative at `alpha_a`.
- `f_b` — function value at `alpha_b`.
- `dir_der_b` — directional derivative at `alpha_b`.

Inputs must be finite real scalars.

## Outputs

- `alpha` — selected maximiser of the cubic model within the interpolation interval.
- `fx` — cubic model value at `alpha`.

## Source

- [Spinach documentation for cubic_interp.m](https://spindynamics.org/wiki/index.php?title=cubic_interp.m)