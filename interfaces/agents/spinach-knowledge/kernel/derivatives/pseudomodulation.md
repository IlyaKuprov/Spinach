# kernel/derivatives/pseudomodulation.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/pseudomodulation.m) · [Spinach Wiki documentation](https://spindynamics.org/wiki/index.php?title=pseudomodulation.m)

## Purpose and signature

`output=pseudomodulation(field,spectrum,mod_amp,mod_order)` applies the requested Fourier-domain pseudomodulation harmonic to one or more uniformly sampled spectra.

## Inputs and output

- `field`: dense floating-point, finite, real `N`-by-1 column with at least three points, uniform spacing, and non-zero sweep width. Either sweep direction is handled.
- `spectrum`: non-empty dense floating-point `N`-by-`M` matrix with finite entries. Rows correspond to field points; each column is an independent spectrum.
- `mod_amp`: finite, non-negative real scalar in the same field units as `field`.
- `mod_order`: scalar harmonic order 0, 1, or 2.
- `output`: same row and column dimensions as `spectrum`; real input spectra are returned real, while complex input may produce a complex result.

The Fourier-domain treatment follows Hyde et al., *Journal of Magnetic Resonance* **96**, 1–13 (1992), Eqs. 5–7, as cited in the MATLAB source.

## Fourier-domain construction

The routine obtains the field step, builds an angular-frequency vector in MATLAB FFT ordering (including the sign of the field step), and forms the dimensionless Bessel argument `mod_amp*ang_freq/2`. It transforms along the row/field dimension, multiplies by the Bessel factor for the selected order, then inverse-transforms along that same dimension:

- Order 0: `ifft(spec_ft.*besselj(0,bessel_arg),npts,1)`.
- Order 1: `2i*ifft(spec_ft.*besselj(1,bessel_arg),npts,1)`.
- Order 2: `2*ifft(spec_ft.*besselj(2,bessel_arg),npts,1)`.

The field spacing supplies the reciprocal-field scale of the frequency axis; the amplitude and field must use matching units so the Bessel argument is dimensionless. The source uses the finite, uniformly spaced axis guard rather than padding or correcting irregular samples.
