# kernel/line_shapes/lorentzfun.m

- Signature: `[real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)`

## Purpose

Normalised Lorentzian function in magnetic resonance notation with a phase distortion. Syntax: [real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

## Parameters / inputs

- offs -peak offset from zero
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension
- phi -phase distortion, radians

## Outputs

- real_part -an array of values, same size as x
- imag_part -an array of values, same size as x

## Implementation structure

- Normalised Lorentzian function in magnetic resonance notation
- with a phase distortion. Syntax:
- [real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)
- offs -peak offset from zero
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension
- phi -phase distortion, radians
- real_part -an array of values, same size as x
- imag_part -an array of values, same size as x
- Check consistency
- Width parameter
