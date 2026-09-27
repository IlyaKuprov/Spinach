# examples/fundamentals/correlation_function_4.m

- Signature: `correlation_function_4()`

## Purpose

Compares a Monte Carlo estimate with the analytical Spinach result for `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>`, using isotropic rotational diffusion at higher tensor rank. The sigma parameter sets the rotational-rate scale; the four indices select the Wigner-function elements being correlated.

## Physical / mathematical content

The test sets `sigma_iso=0.2`, rank `L=4`, and indices `k=-1, m=2, p=-1, q=2` (converted from `[-L,L]` indexing to MATLAB array indices). The analytical model uses Redfield relaxation and correlation time `1/(3*sigma_iso^2)`.

## Numerical / algorithmic content

It generates `1e6` rotations and estimates the correlation over `nlags=100`. Direction-cosine matrices are converted to Wigner functions, and normalized cross-correlation of the selected elements is compared with the exponential sum from Spinach's `corrfun` calculation.

## Implementation structure

The Monte Carlo calculation stores the rotation trajectory and evaluates Wigner matrices in a `parfor` loop. A one-spin dummy system supplies the analytical curve; the example plots both results. The source estimates a run time of minutes.
