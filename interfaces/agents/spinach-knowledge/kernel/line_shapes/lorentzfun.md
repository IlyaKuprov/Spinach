# kernel/line_shapes/lorentzfun.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/line_shapes/lorentzfun.m`
- Signature: `[real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)`
- Total lines: 69

## Purpose

Normalised Lorentzian function in magnetic resonance notation with a phase distortion. Syntax: [real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(offs,ampl,fwhm,x,phi)`.
- Lines 34-35: Width parameter; implemented by `gam=fwhm/2`.
- Lines 37-39: Calculate output; implemented by `real_part=((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*cos(phi)- ((x-offs)/gam).*((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*sin(phi)`.

### Key state/data transformations

- Lines 35: computes `gam` using `gam=fwhm/2`.
- Lines 38-39: computes `real_part` using `real_part=((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*cos(phi)- ((x-offs)/gam).*((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*sin(phi)`.
- Lines 40-41: computes `imag_part` using `imag_part=((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*sin(phi)+ ((x-offs)/gam).*((ampl/(2*pi*gam))./(1+((x-offs)/gam).^2)).*cos(phi)`.

### Local helper functions

- Line 46: `grumble()` — `function grumble(offs,ampl,fwhm,x,phi)`.
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))`.
  - Representative operation: `error('x must be an array of real numbers.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
