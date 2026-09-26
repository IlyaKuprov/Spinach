# kernel/line_shapes/gausscon.m

- Signature: `y=gausscon(offs,ampl,fwhm,x)`

## Purpose

Normalised Gaussian function in magnetic resonance notation and its convolution with a triangular function. Syntax: y=gausscon(offs,ampl,fwhm,x)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

## Parameters / inputs

- offs -peak offset from zero -when this is a scalar,
- a Gaussian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension

## Outputs

- y -an array of values, same size as x

## Implementation structure

- Normalised Gaussian function in magnetic resonance notation and
- its convolution with a triangular function. Syntax:
- y=gausscon(offs,ampl,fwhm,x)
- offs -peak offset from zero -when this is a scalar,
- a Gaussian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension
- y -an array of values, same size as x
- Check consistency
