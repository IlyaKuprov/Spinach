# kernel/line_shapes/gaussfun.m

- Signature: `y=gaussfun(x,fwhm)`

## Purpose

Evaluates a Gaussian line shape centered at zero and normalized to unit area.

## Physical / mathematical content

The width is specified by the full width at half maximum, `fwhm`.

## Numerical / algorithmic content

The routine converts `fwhm` to the Gaussian standard deviation and evaluates the normalized Gaussian elementwise at `x`.

## Parameters / inputs

- `x` - real numeric array of any dimension.
- `fwhm` - positive real scalar full width at half maximum.

## Outputs

- `y` - Gaussian values at the points in `x`, with the same array dimensions.
