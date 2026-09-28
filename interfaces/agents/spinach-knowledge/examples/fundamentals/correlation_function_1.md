# examples/fundamentals/correlation_function_1.m

- Signature: `correlation_function_1()`

## Purpose

Compares a Monte Carlo estimate with the analytical Spinach result for the rotational correlation function `G(L,k,m,p,q)=<D{L}(k,m)*D{L}(p,q)'>`. The example tests isotropic rotational diffusion; `sigma_iso` sets the rotational-rate scale, and the indices select the correlated Wigner-function elements.

## Physical / mathematical content

The test uses `sigma_iso=0.1`, rank `L=2`, and indices `k=-1, m=2, p=-1, q=2` (converted from `[-L,L]` indexing to MATLAB array indices). The analytical model is represented with Redfield relaxation and correlation time `1/(3*sigma_iso^2)`.

## Numerical / algorithmic content

A Monte Carlo trajectory of `1e6` rotations is generated from Gaussian increments, with `nlags=300`. Each direction-cosine matrix is converted to Euler angles and then to Wigner functions; normalized cross-correlation is compared with the exponential sum returned by Spinach's `corrfun` calculation.

## Implementation structure

The script constructs and stores the rotation-matrix trajectory, computes the selected Wigner-element cross-correlation, builds a one-spin dummy system for the analytical calculation, and plots Monte Carlo points against the Spinach curve. The source estimates a run time of minutes.
