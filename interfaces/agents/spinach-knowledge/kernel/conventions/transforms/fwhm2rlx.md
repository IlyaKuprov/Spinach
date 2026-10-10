# kernel/conventions/transforms/fwhm2rlx.m

## Purpose

Approximates the transverse relaxation rate from an NMR signal's full width at half-maximum (FWHM), assuming a Lorentzian line shape.

## Signature

`r2rate=fwhm2rlx(fwhm)`

## Conversion

`r2rate = pi * fwhm`

Both the input FWHM and output approximate R2 rate are in Hz.

## Input and output

- `fwhm`: real numeric array. The function errors if any element is less than or equal to zero. The validation does not explicitly test finiteness.
- `r2rate`: numeric array with the same dimensions as `fwhm`, scaled elementwise by pi.

FWHM alone is not a reliable measure of transverse relaxation; the source says to treat the result as an upper bound. The Lorentzian line shape is assumed.

## References

- MATLAB source: [kernel/conventions/transforms/fwhm2rlx.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/fwhm2rlx.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=fwhm2rlx.m)
