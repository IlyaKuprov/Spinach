# kernel/utilities/cheb_coeff.m

- Signature: `c=cheb_coeff(f,a,b,n)`

## Purpose

Discrete cosine transform algorithm for Chebyshev expansion coefficients of the user-specified scalar function. Syntax: c=cheb_coeff(f,a,b,n)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

## Parameters / inputs

- f -function handle, must be vectorised
- a -left edge of the expansion interval
- b -right edge of the expansion interval
- n -number of Chebyshev polynomials in
- the expansion

## Outputs

- c -a vector of expansion coefficients

## Implementation structure

- Discrete cosine transform algorithm for Chebyshev expansion
- coefficients of the user-specified scalar function. Syntax:
- c=cheb_coeff(f,a,b,n)
- f -function handle, must be vectorised
- a -left edge of the expansion interval
- b -right edge of the expansion interval
- n -number of Chebyshev polynomials in
- the expansion
- c -a vector of expansion coefficients
- Check consistency
- [-1,+1] query points
- Scaled query points
