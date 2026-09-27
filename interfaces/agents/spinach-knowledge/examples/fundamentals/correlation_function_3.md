# examples/fundamentals/correlation_function_3.m

- Signature: `correlation_function_3()`

## Purpose

Compares a Monte Carlo estimate with the analytical Spinach result for `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>` under rhombic rotational diffusion. The sigma parameters set the three rotational-rate scales; the indices identify the Wigner-function elements being correlated.

## Physical / mathematical content

The example uses `sigma_x=0.1`, `sigma_y=0.2`, `sigma_z=0.3`, rank `L=2`, and indices `k=-1, m=2, p=-1, q=2` (converted from `[-L,L]` indexing to MATLAB array indices). The analytical model uses Redfield relaxation and correlation times `1./(3*[sigma_x sigma_y sigma_z].^2)`.

## Numerical / algorithmic content

A Monte Carlo trajectory of `1e6` rotations is generated, and the estimate uses `nlags=300`. Independent Gaussian increments are scaled by the three sigma values; the resulting direction-cosine matrices are converted to Wigner functions and the selected elements are cross-correlated. This is compared to the exponential sum returned by Spinach's `corrfun` calculation.

## Implementation structure

The code stores the rotation trajectory, computes Wigner matrices with `parfor`, builds the dummy system for the analytical calculation, and plots both curves. The source estimates a run time of minutes.
