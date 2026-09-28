# kernel/line_shapes/dhofun.m

- Signature: `y=dhofun(x,nat_freq,fwhm)`

## Purpose

Evaluates the normalized damped-harmonic-oscillator response at positive frequencies. It peaks at the undamped natural frequency and goes to zero quadratically as frequency approaches zero.

## Physical / mathematical content

The response models a damped oscillator band. In the weak-damping limit it approaches a Lorentzian; `fwhm` is the full width at half maximum for all damping values.

## Numerical / algorithmic content

Non-positive entries of `x` return zero. The formula is evaluated using frequencies scaled by `nat_freq` to avoid intermediate overflow in single precision.

## Parameters / inputs

- `x` - real numeric array of any dimension, in frequency units.
- `nat_freq` - finite positive real scalar, the natural frequency of the undamped oscillator, in the same units as `x`.
- `fwhm` - finite positive real scalar, the oscillator damping rate and full width at half maximum, in the same units as `x`.

## Outputs

- `y` - response values with the same size and type as `x`, normalized to unit integral over positive frequencies.
