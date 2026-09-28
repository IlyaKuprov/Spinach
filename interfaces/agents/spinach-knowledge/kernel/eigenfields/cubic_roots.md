# kernel/eigenfields/cubic_roots.m

- Signature: `root_list=cubic_roots(poly_coeffs,root_tol)`
- Source URL: https://spindynamics.org/wiki/index.php?title=cubic_roots.m

## Purpose

Find real roots of a cubic polynomial in the unit interval.

## Parameters / inputs

- `poly_coeffs`: four finite real coefficients `[a b c d]` of `a*x^3+b*x^2+c*x+d`.
- `root_tol`: finite positive real scalar used as a root-filtering tolerance.

## Output

- `root_list`: sorted row vector of accepted real roots in `[0,1]`; empty if no roots are returned.

## Numerical / algorithmic content

1. Validate both inputs, reshape the coefficients to a row vector, and divide them by their maximum absolute value. Return an empty array if that maximum is zero.
2. Drop leading coefficients whose absolute values do not exceed `root_tol`. Return an empty array if none remain.
3. Compute roots of the remaining polynomial. Retain roots only when `abs(imag(root)) < root_tol`, then take their real parts.
4. Retain real parts in `[-root_tol,1+root_tol]`, clamp them to `[0,1]`, and sort them into a row vector.
5. Merge adjacent sorted roots whose difference does not exceed `root_tol`.