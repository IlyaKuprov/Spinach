# kernel/eigenfields/cubic_roots.m

- Direct source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/eigenfields/cubic_roots.m
- Spinach Wiki: https://spindynamics.org/wiki/index.php?title=cubic_roots.m

## Signature

`root_list=cubic_roots(poly_coeffs,root_tol)`

## Purpose

Find the real roots of a polynomial of degree at most three in the unit interval. This helper is used by eigenfield calculations to solve cubic interpolants on a normalised field-interval coordinate.

## Inputs

- `poly_coeffs`: numeric real array with four finite elements, interpreted as `[a b c d]` for `a*x^3+b*x^2+c*x+d`. The input is reshaped to a row.
- `root_tol`: finite positive real numeric scalar used both to discard small leading coefficients and to filter roots.

## Output

- `root_list`: sorted row vector of accepted real roots in `[0,1]`; empty when no root is accepted.

## Algorithm and edge cases

The coefficients are divided by their maximum absolute value, making the computation insensitive to a common nonzero coefficient scale. All-zero coefficients return an empty result. Leading coefficients with absolute value at most `root_tol` are dropped, so a lower-degree polynomial is handled; if no coefficient remains, the result is empty. MATLAB's polynomial root calculation is then filtered by `abs(imag(root)) < root_tol`. Real parts within `[-root_tol,1+root_tol]` are clamped to `[0,1]`, sorted, and adjacent roots separated by at most `root_tol` are merged.

The polynomial coordinate is dimensionless and the output is in the unit interval. The source contains no worked numeric example.
