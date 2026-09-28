# examples/fundamentals/correlation_function_2.m

- Signature: `correlation_function_2()`

## Purpose

Compares a Monte Carlo estimate with the analytical Spinach result for `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>`. This example tests axial rotational diffusion; the sigma parameters specify the rates along the axial and equivalent directions, while the four indices select the Wigner-function elements.

## Physical / mathematical content

The parameters are `sigma_ax=0.1`, `sigma_eq=0.2`, rank `L=2`, and indices `k=-1, m=2, p=-1, q=2` (converted from `[-L,L]` indexing to MATLAB array indices). The analytical model uses Redfield relaxation with correlation times `1./(3*[sigma_ax sigma_eq].^2)`.

## Numerical / algorithmic content

The Monte Carlo trajectory contains `1e6` rotations and the estimate uses `nlags=300`. Gaussian increments are scaled separately for axial and equivalent rotations; direction-cosine matrices are converted to Wigner functions, whose selected elements are cross-correlated. The normalized estimate is compared to the exponential sum from Spinach's `corrfun` calculation.

## Implementation structure

The source stores the trajectory, evaluates Wigner matrices in a `parfor` loop, then computes both the Monte Carlo correlation and the analytical Spinach curve for plotting. It labels the source's estimated run time as minutes.
