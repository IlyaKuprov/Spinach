# kernel/line_shapes/dhofun.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/line_shapes/dhofun.m`
- Signature: `y=dhofun(x,nat_freq,fwhm)`
- Total lines: 65

## Purpose

Normalised damped harmonic oscillator response function in mag- netic resonance notation. This is the standard shape of a phonon band: a resonance at the natural frequency of the oscillator, br- oadened by damping, and vanishing quadratically at zero frequen- cy, which a Lorentzian centred at the same place does not do. In the weak damping limit the function tends to a Lorentzian of the same width. Evaluation uses fr

## Physical / mathematical content

- Line-shape utilities. These files compute, transform, or fit spectral line shapes, connecting simulated transition frequencies and relaxation widths to observable spectra.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(x,nat_freq,fwhm)`.
- Lines 41-42: The response function only lives at positive frequencies; implemented by `y=zeros(size(x),'like',x); pos_args=(x>0); arg=x(pos_args)`.
- Lines 44-45: Compute the response in reduced frequency units; implemented by `rel_freq=arg/nat_freq; rel_fwhm=fwhm/nat_freq`.

### Key state/data transformations

- Lines 42: computes `y` using `y=zeros(size(x),'like',x); pos_args=(x>0); arg=x(pos_args)`.
- Lines 45: computes `rel_freq` using `rel_freq=arg/nat_freq; rel_fwhm=fwhm/nat_freq`.
- Lines 46-47: computes `y(pos_args)` using `y(pos_args)=(2*rel_fwhm/(pi*nat_freq))*(rel_freq.^2)./ ((rel_freq.^2-1).^2+(rel_fwhm*rel_freq).^2)`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(x,nat_freq,fwhm)`.
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))`.
  - Representative operation: `error('x must be an array of real numbers.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`.
