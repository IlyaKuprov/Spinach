# kernel/line_shapes/gaussfun.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/line_shapes/gaussfun.m`
- Signature: `y=gaussfun(x,fwhm)`
- Total lines: 49

## Purpose

Normalized Gaussian function in magnetic resonance notation. Syntax: y=gaussfun(x,fwhm)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 22-23: Check consistency; implemented by `grumble(x,fwhm)`.
- Lines 25-26: Compute standard deviation; implemented by `sigma=fwhm/(2*sqrt(2*log(2)))`.
- Lines 28-29: Compute the Gaussian; implemented by `y=(1/(sigma*sqrt(2*pi)))*exp(-(x.^2)/(2*sigma^2))`.

### Key state/data transformations

- Lines 26: computes `sigma` using `sigma=fwhm/(2*sqrt(2*log(2)))`.
- Lines 29: computes `y` using `y=(1/(sigma*sqrt(2*pi)))*exp(-(x.^2)/(2*sigma^2))`.

### Local helper functions

- Line 34: `grumble()` — `function grumble(x,fwhm)`. Fifty years ago the back streets of Leningrad have taught me one lesson: when a fight is un-
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))`.
  - Representative operation: `error('x must be an array of real numbers.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
