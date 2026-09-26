# kernel/line_shapes/dhofun.m

- Signature: `y=dhofun(x,nat_freq,fwhm)`

## Purpose

Normalised damped harmonic oscillator response function in mag- netic resonance notation. This is the standard shape of a phonon band: a resonance at the natural frequency of the oscillator, br- oadened by damping, and vanishing quadratically at zero frequen- cy, which a Lorentzian centred at the same place does not do. In the weak damping limit the function tends to a Lorentzian of the same width. Evaluation uses fr

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

## Parameters / inputs

- x -argument values, a real array of any dimension;
- the function is zero at non-positive arguments
- nat_freq -natural frequency of the undamped oscillator, a
- positive real scalar, in the same units as x; the
- maximum of the function sits exactly here
- fwhm -damping rate of the oscillator, a positive real
- scalar, in the same units as x; for this respon-
- se function it is exactly the full width at half-
- maximum, at any damping

## Outputs

- y -function values at the points specified in x, an
- array of the same size and type as x, normalised
- to unit integral over positive arguments

## Implementation structure

- Normalised damped harmonic oscillator response function in mag-
- netic resonance notation. This is the standard shape of a phonon
- band: a resonance at the natural frequency of the oscillator, br-
- oadened by damping, and vanishing quadratically at zero frequen-
- cy, which a Lorentzian centred at the same place does not do. In
- the weak damping limit the function tends to a Lorentzian of the
- same width. Evaluation uses frequencies scaled by nat_freq to
- avoid intermediate overflow in single precision. Syntax:
- y=dhofun(x,nat_freq,fwhm)
- x -argument values, a real array of any dimension;
- the function is zero at non-positive arguments
- nat_freq -natural frequency of the undamped oscillator, a
