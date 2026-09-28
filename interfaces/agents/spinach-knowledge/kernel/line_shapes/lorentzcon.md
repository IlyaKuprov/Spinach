# kernel/line_shapes/lorentzcon.m

- Signature: `y=lorentzcon(offs,ampl,fwhm,x)`

## Purpose

Evaluates a Lorentzian line shape or its convolution with a boxcar or triangular distribution.

## Physical / mathematical content

A scalar `offs` gives a Lorentzian centered at that offset. Two offsets specify a boxcar distribution; three offsets specify the vertices of a triangular distribution. `ampl` scales the resulting values.

## Numerical / algorithmic content

The routine evaluates the Lorentzian or its analytic convolution, sorting the offsets and handling repeated triangle vertices. Integer-valued `x` is converted to double precision before evaluation.

## Parameters / inputs

- `offs` - finite real vector with one, two, or three elements.
- `ampl` - finite real scalar amplitude multiplier.
- `fwhm` - finite positive real scalar full width at half maximum.
- `x` - array of finite real argument values.

## Outputs

- `y` - array of values with the same size as `x`.
