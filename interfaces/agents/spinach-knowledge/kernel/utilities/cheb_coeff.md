# kernel/utilities/cheb_coeff.m

- Signature: `c=cheb_coeff(f,a,b,n)`

## Purpose

Computes the coefficients of a Chebyshev expansion of the user-specified scalar function over the interval `[a,b]` using a discrete cosine transform.

## Physical / mathematical content

The function samples `f` at Chebyshev query points mapped from `[-1,1]` to `[a,b]`, then returns the corresponding expansion coefficients.

## Numerical / algorithmic content

After checking the inputs, it computes `dct(f(x))/sqrt(n)` and scales coefficients 2 through `n` by `sqrt(2)`.

## Parameters / inputs

- `f` — vectorised function handle
- `a` — left edge of the expansion interval
- `b` — right edge of the expansion interval
- `n` — number of Chebyshev polynomials in the expansion

## Outputs

- `c` — vector of expansion coefficients
