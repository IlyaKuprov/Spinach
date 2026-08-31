# kernel/utilities/cheb_coeff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/cheb_coeff.m`
- Signature: `c=cheb_coeff(f,a,b,n)`
- Total lines: 61

## Purpose

Discrete cosine transform algorithm for Chebyshev expansion coefficients of the user-specified scalar function. Syntax: c=cheb_coeff(f,a,b,n)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(f,a,b,n)`.
- Lines 30-31: [-1,+1] query points; implemented by `x=cos(((1:n)*2-1)*pi/(2*n))`.
- Lines 33-34: Scaled query points; implemented by `x=0.5*(a+x*(b-a)+b)`.
- Lines 36-37: Expansion coefficients; implemented by `c=dct(f(x))/sqrt(n)`.

### Key state/data transformations

- Lines 31: computes `x` using `x=cos(((1:n)*2-1)*pi/(2*n))`.
- Lines 37: computes `c` using `c=dct(f(x))/sqrt(n)`.
- Lines 38: computes `c(2:n)` using `c(2:n)=c(2:n)*sqrt(2)`.

### Local helper functions

- Line 43: `grumble()` — `function grumble(f,a,b,n)`.
  - Representative operation: `if ~isa(f,'function_handle')`.
  - Representative operation: `error('f must be a function handle.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dct()`, `isscalar()`.
