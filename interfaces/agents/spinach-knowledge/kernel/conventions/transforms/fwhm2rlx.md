# kernel/conventions/transforms/fwhm2rlx.m

- Signature: `r2rate=fwhm2rlx(fwhm)`

## Purpose

Approximates the transverse relaxation rate `R2` from an NMR signal's full width at half-maximum, assuming a Lorentzian line shape. The result should be treated as an upper bound: FWHM is not a reliable measure of transverse relaxation by itself.

## Parameters / inputs

- `fwhm`: positive real value or array of linewidths in Hz.

## Output

- `r2rate`: approximate `R2` rate in Hz, calculated as `pi * fwhm`.

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=fwhm2rlx.m)
