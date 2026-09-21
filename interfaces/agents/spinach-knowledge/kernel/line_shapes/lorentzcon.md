# kernel/line_shapes/lorentzcon.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/line_shapes/lorentzcon.m`
- Signature: `y=lorentzcon(offs,ampl,fwhm,x)`
- Total lines: 106

## Purpose

Normalised Lorentzian function in magnetic resonance notation and its convolution with a triangular function. Syntax: y=lorentzcon(offs,ampl,fwhm,x)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- offs -peak offset from zero -when this is a scalar,
- a Lorentzian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension

## Outputs

- y -an array of values, same size as x

## Implementation structure

- Normalised Lorentzian function in magnetic resonance notation and
- its convolution with a triangular function. Syntax:
- y=lorentzcon(offs,ampl,fwhm,x)
- offs -peak offset from zero -when this is a scalar,
- a Lorentzian is returned; when this is a vector
- with three elements, a convolution with a tri-
- angular function is returned.
- ampl -amplitude multiplier, scalar
- fwhm -full width at half-maximum, scalar
- x -argument, array of any dimension
- y -an array of values, same size as x
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `offs()`, `isscalar()`, `elseif()`, `atan2()`, `ismember()`, `any()`.
