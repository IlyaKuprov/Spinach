# kernel/utilities/cheb_coeff.m

## Purpose

Computes the Chebyshev expansion coefficients of a user-specified scalar function using a discrete cosine transform (DCT) algorithm.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cheb_coeff.m>

## Behaviour

- Syntax: `c=cheb_coeff(f,a,b,n)`.
- Validates inputs via an internal `grumble` subfunction, which errors if `f` is not a function handle, if `a` and `b` are not real scalars with `a < b`, or if `n` is not a positive real integer.
- Generates `n` Chebyshev–Gauss query points on `[-1,+1]` as `x=cos(((1:n)*2-1)*pi/(2*n))`.
- Scales the query points to the interval `[a,b]` via `x=0.5*(a+x*(b-a)+b)`.
- Evaluates the function at the scaled points and applies `dct`, dividing by `sqrt(n)`; coefficients `c(2:n)` are then multiplied by `sqrt(2)`.

## Inputs and outputs

Inputs:

- `f` — function handle, must be vectorised.
- `a` — left edge of the expansion interval.
- `b` — right edge of the expansion interval.
- `n` — number of Chebyshev polynomials in the expansion.

Output:

- `c` — a vector of expansion coefficients.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=cheb_coeff.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/cheb_coeff.m>
