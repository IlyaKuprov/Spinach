# examples/fundamentals/correlation_function_5.m

- Signature: `correlation_function_5()`

## Purpose

Uses Monte Carlo rotational diffusion to estimate `G(k,m,p,q)=<R(k,m)*R(p,q)>`, where `R` is a three-dimensional Cartesian rotation matrix. Unlike the preceding Wigner-function examples, this script plots only the Monte Carlo estimate.

## Physical / mathematical content

The isotropic rate parameter is `sigma_iso=0.2`; the selected matrix elements are `k=2, m=3, p=2, q=3`. The correlation is scaled by `1/3`, as implemented in the source.

## Numerical / algorithmic content

The script propagates `1e6` rotations from Gaussian angular increments and computes a normalized cross-correlation with `nlags=300`. It plots the real Monte Carlo correlation against lag; the source estimates a run time of minutes.

## Implementation structure

Starting from the identity matrix, each step right-multiplies the accumulated rotation by the matrix exponential of the increment generator scaled by `sigma_iso`. The chosen Cartesian matrix elements are passed to `xcorr`; the shifted result is scaled by `1/3` before plotting.
