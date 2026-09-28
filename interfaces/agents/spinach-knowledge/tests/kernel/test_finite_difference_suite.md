# tests/kernel/test_finite_difference_suite.m

- Signature: `result=test_finite_difference_suite()`

## Purpose

Checks finite-difference and spectral-differentiation helpers against exact or analytically known cases.

## Physical / mathematical content

The suite tests numerical derivative identities, including finite-difference and Fourier differentiation, periodic operators, and directional derivatives of a commuting matrix exponential. It is a helper regression suite rather than acquired-spectrum processing.

## Numerical / algorithmic content

- Checks interpolation and first- and second-derivative finite-difference weights, finite-difference matrices, Fourier differentiation, Laplacians, FFT differentiation kernels, pseudomodulation, and matrix-exponential directional derivatives.
- Three-point centred weights at zero reproduce the exact interpolation and derivative coefficients; a five-point wall matrix differentiates quadratics on a unit grid; a seven-point, cubic Savitzky-Golay derivative recovers a cubic exactly.
- Periodic finite-difference and Laplacian matrices annihilate constants. The suite also checks Savitzky-Golay window constraints and pseudomodulation axis conventions.

## Outputs

`result` is the regression-test result with explanatory messages.

## Implementation structure

The test initializes a regression result and evaluates each helper against an exact or closed-form reference, ending with zeroth- and first-directional-derivative checks for a commuting matrix exponential.
