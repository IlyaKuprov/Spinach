# kernel/line_shapes/gausscon.m

- Signature: `y=gausscon(offs,ampl,fwhm,x)`

## Purpose

Evaluates a Gaussian line shape or its convolution with a triangular distribution.

## Physical / mathematical content

A scalar `offs` gives a Gaussian centered at that offset. Three offsets specify the vertices of a triangular distribution whose convolution with the Gaussian is evaluated. `ampl` scales the resulting values.

## Numerical / algorithmic content

The Gaussian standard deviation is obtained from `fwhm`. For three offsets, the routine sorts them and evaluates the convolution using Gaussian values and error-function integrals; repeated vertices are handled as limiting cases.

## Parameters / inputs

- `offs` - finite real scalar or three-element vector of offsets.
- `ampl` - finite real scalar amplitude multiplier.
- `fwhm` - finite positive real full width at half maximum.
- `x` - array of finite real argument values.

## Outputs

- `y` - array of values with the same size as `x`.
