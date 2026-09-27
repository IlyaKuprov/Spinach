# etc/data_processing/autophase.m

- Signature: `[spec,cheb_coeffs]=autophase(spec,guess)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=autophase.m)

## Purpose

Correct the phase of a one-dimensional complex NMR spectrum by fitting a smooth phase profile represented by Chebyshev polynomials.

## Method

The spectral window is mapped to `[-1,1]`. `autophase` optimises the Chebyshev coefficients with `fminunc`, choosing phases that move the spectrum’s fourth-norm signal from the imaginary component toward the real component. The fitted phase multipliers are applied to the spectrum; the initial standard-deviation scaling is then undone.

## Inputs

- `spec` — finite, non-constant numeric vector containing the complex spectrum.
- `guess` — finite real row vector of at least two initial Chebyshev coefficients, in radians. `[phi 0 0]` is a suggested initial value, where `phi` is the zero-order phase guess.

## Outputs

- `spec` — phase-corrected spectrum, returned as a column vector.
- `cheb_coeffs` — fitted Chebyshev coefficients describing the phase profile across the spectral window.
