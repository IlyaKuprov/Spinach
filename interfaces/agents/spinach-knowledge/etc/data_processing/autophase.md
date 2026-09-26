# etc/data_processing/autophase.m

- Signature: `[spec,cheb_coeffs]=autophase(spec,guess)`

## Purpose

Chebyshev phase corrector for 1D NMR spectra. Views the phase profile across the spectral window as a slowly va- rying function and approximates it with a linear combi- nation of low-order Chebyshev polynomials. Syntax: [spec,cheb_coeffs]=autophase(spec,guess)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- spec -1D NMR spectrum, a complex vector
- guess -initial guess for the Chebyshev polynomi-
- al coefficients, radians. [phi 0 0] is a
- good start, where phi is the zero-order
- phase correction guess.

## Outputs

- spec -phased NMR spectrum, a column vector
- coeffs -Chebyshev polynomial coefficients of the
- phase profile across the spectrum with
- the window treated as a [-1,1] interval

## Implementation structure

- Chebyshev phase corrector for 1D NMR spectra. Views the
- phase profile across the spectral window as a slowly va-
- rying function and approximates it with a linear combi-
- nation of low-order Chebyshev polynomials. Syntax:
- [spec,cheb_coeffs]=autophase(spec,guess)
- spec -1D NMR spectrum, a complex vector
- guess -initial guess for the Chebyshev polynomi-
- al coefficients, radians. [phi 0 0] is a
- good start, where phi is the zero-order
- phase correction guess.
- spec -phased NMR spectrum, a column vector
- coeffs -Chebyshev polynomial coefficients of the
