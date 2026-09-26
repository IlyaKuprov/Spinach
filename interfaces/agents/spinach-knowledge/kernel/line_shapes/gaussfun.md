# kernel/line_shapes/gaussfun.m

- Signature: `y=gaussfun(x,fwhm)`

## Purpose

Normalized Gaussian function in magnetic resonance notation. Syntax: y=gaussfun(x,fwhm)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

## Parameters / inputs

- x -argument values, a real array of any dimension
- fwhm -full width at half-maximum

## Outputs

- y -function values at the points specified in x

## Implementation structure

- Normalized Gaussian function in magnetic resonance
- notation. Syntax:
- y=gaussfun(x,fwhm)
- x -argument values, a real array of any dimension
- fwhm -full width at half-maximum
- y -function values at the points specified in x
- Check consistency
- Compute standard deviation
- Compute the Gaussian
- Consistency enforcement
- Fifty years ago the back streets of Leningrad
- have taught me one lesson: when a fight is un-
