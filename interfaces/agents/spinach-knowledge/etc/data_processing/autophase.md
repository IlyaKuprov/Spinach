# etc/data_processing/autophase.m

- MATLAB implementation: [etc/data_processing/autophase.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/data_processing/autophase.m)

## Purpose and call

Correct the phase of a one-dimensional complex NMR spectrum by fitting a phase profile across the spectral window with low-order Chebyshev polynomials.

`[spec,cheb_coeffs]=autophase(spec,guess)`

`spec` must be a finite numeric vector with non-zero standard deviation. The initial `guess` must be a finite real row vector of at least two coefficients, in radians; `[phi 0 0]` is the suggested starting form, with `phi` the zero-order phase guess.

## Numerical mechanism

The routine divides the spectrum by `std(spec)` before optimisation, then uses `fminunc` with central finite differences, no iteration/evaluation cap, and function/optimality/step tolerances of `1e-12`. For sample positions mapped linearly across `[-1,1]`, it builds Chebyshev polynomials by recurrence and applies the pointwise phase multiplier `exp(1i*phis*cheb)`. The objective minimised is `norm(imag(spec),4)-norm(real(spec),4)`, moving fourth-norm signal from the imaginary to the real component. It applies the fitted correction, restores the original standard-deviation scale, and returns the coefficients.

## Outputs and limitation

`spec` is returned as the phased spectrum in a column vector; `cheb_coeffs` are the fitted phase-profile coefficients for the window mapped to `[-1,1]`. The fit is an unconstrained `fminunc` optimisation initialised by `guess`, so the starting coefficients matter; the function does not impose phase bounds.

Source and credit: [Spinach Wiki: autophase.m](https://spindynamics.org/wiki/index.php?title=autophase.m). Ilya Kuprov (contact detail is in the source file).