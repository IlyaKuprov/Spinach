# etc/data_processing/autophase.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/data_processing/autophase.m`
- Signature: `[spec,cheb_coeffs]=autophase(spec,guess)`
- Total lines: 98

## Purpose

Chebyshev phase corrector for 1D NMR spectra. Views the phase profile across the spectral window as a slowly va- rying function and approximates it with a linear combi- nation of low-order Chebyshev polynomials. Syntax: [spec,cheb_coeffs]=autophase(spec,guess)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `objective()`, `apply_phases()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `std()`, `optimoptions()`, `fminunc()`, `objective()`, `apply_phases()`, `cheb()`, `phi_mults()`, `spec()`, `isvector()`, `any()`, `isrow()`, `guess()`.
