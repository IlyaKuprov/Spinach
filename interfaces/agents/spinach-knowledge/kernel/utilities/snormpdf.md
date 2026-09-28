# kernel/utilities/snormpdf.m

- Signature: `p=snormpdf(x,mu,sigma,alpha)`

## Purpose

Compute Azzalini's skew normal probability density.

## Physical / mathematical content

- `p=2*normpdf(x,mu,sigma).*normcdf(alpha*x,alpha*mu,sigma)` (Equation 2 in http://www.jstor.org/stable/4615982).
- `mu` is the expectation value of the normal distribution.

## Numerical / algorithmic content

- Evaluates the density elementwise for `x`.

## Parameters / inputs

- `x` — an array of real numbers.
- `mu` — expectation value of the normal distribution; a real scalar.
- `sigma` — standard deviation of the normal distribution; a positive real scalar.
- `alpha` — skew factor; a real scalar.

## Outputs

- `p` — an array of probability densities with the same shape as `x`.

## Implementation structure

- Checks input consistency, then evaluates the density formula.
- Source: https://spindynamics.org/wiki/index.php?title=snormpdf.m
- Contact: ilya.kuprov@weizmann.ac.il

- Equation 2 in http://www.jstor.org/stable/4615982

<https://spindynamics.org/wiki/index.php?title=snormpdf.m>