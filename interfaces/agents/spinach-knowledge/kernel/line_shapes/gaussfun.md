# kernel/line_shapes/gaussfun.m

- MATLAB source: [kernel/line_shapes/gaussfun.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/line_shapes/gaussfun.m)
- Existing Wiki: [gaussfun.m](https://spindynamics.org/wiki/index.php?title=gaussfun.m)
- Signature: `y=gaussfun(x,fwhm)`

## Meaning and equation

This evaluates a zero-centred Gaussian normalised to unit area. The source sets `sigma=fwhm/(2*sqrt(2*log(2)))` and evaluates `y=exp(-x^2/(2*sigma^2))/(sigma*sqrt(2*pi))` elementwise. Thus `fwhm` is the full width at half maximum, and the integral over the real line is one.

## Inputs and units

- `x` - real numeric array of any dimension. The guard checks numeric and real input but does not explicitly require finite entries.
- `fwhm` - positive real numeric scalar. The source checks the scalar count and positivity; it does not explicitly reject non-finite values.

`x` and `fwhm` use the same coordinate units. The function performs no conversion or distinction between Hz and angular frequency; use one convention consistently. The line-shape values have reciprocal-coordinate units.

## Output

- `y` - Gaussian values with the same size as `x`.
