# kernel/derivatives/pseudomodulation.m

- Signature: `output=pseudomodulation(field,spectrum,mod_amp,mod_order)`

## Purpose

Pseudomodulation of uniformly sampled spectra using the Hyde et al. Fourier-domain algorithm.

## Parameters / inputs

- `field` — N-by-1 real, ordered, uniformly spaced field axis with at least three points and a non-zero sweep width.
- `spectrum` — N-by-M spectrum matrix; rows are field samples and columns are independent spectra.
- `mod_amp` — non-negative, finite real scalar modulation amplitude in field units.
- `mod_order` — modulation harmonic order: 0, 1, or 2.

## Output

- `output` — N-by-M pseudomodulated spectrum matrix. For real input spectra, round-off imaginary parts are removed.

## Numerical / algorithmic content

The implementation follows Eqs. 5–7 of Hyde et al., J. Magn. Reson. 96, 1–13 (1992). It constructs an angular-frequency axis in MATLAB FFT ordering from the field spacing, then sets the Bessel-function argument to `mod_amp*ang_freq/2`. It Fourier-transforms each spectrum along the field axis, multiplies by the Bessel function for the requested harmonic, and inverse-transforms:

- Order 0: `ifft(spec_ft.*besselj(0,bessel_arg),npts,1)`.
- Order 1: `2i*ifft(spec_ft.*besselj(1,bessel_arg),npts,1)`.
- Order 2: `2*ifft(spec_ft.*besselj(2,bessel_arg),npts,1)`.

After phase-sensitive detection, the time-dependent prefactors are set to unity, leaving amplitude factors `2i` for the first harmonic and `2` for the second harmonic.

[Source documentation](https://spindynamics.org/wiki/index.php?title=pseudomodulation.m)