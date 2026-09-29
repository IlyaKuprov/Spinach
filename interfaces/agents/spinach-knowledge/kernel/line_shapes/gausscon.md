# kernel/line_shapes/gausscon.m

- MATLAB source: [kernel/line_shapes/gausscon.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/line_shapes/gausscon.m)
- Existing Wiki: [gausscon.m](https://spindynamics.org/wiki/index.php?title=gausscon.m)
- Signature: `y=gausscon(offs,ampl,fwhm,x)`

## Meaning and equation

This convolves a unit-area Gaussian with either a point offset or a normalised triangular distribution. The Gaussian kernel is `g(u)=exp(-u^2/(2*sigma^2))/(sigma*sqrt(2*pi))`, with `sigma=fwhm/(2*sqrt(2*log(2)))`.

- One offset `q`: `y(x)=ampl*g(x-q)`.
- Three distinct sorted offsets `a<b<c`: use the triangular density `T(t)=2*(t-a)/((b-a)*(c-a))` for `a<=t<=b`, `T(t)=2*(c-t)/((c-b)*(c-a))` for `b<t<=c`, and zero elsewhere. Then `y(x)=ampl*integral(T(t)*g(x-t),t=a..c)`.

The source sorts the three offsets. Repeated or numerically coalescent vertices are handled as the corresponding limiting shape: all coalescent gives the Gaussian at their mean; a repeated lower pair or upper pair gives a right-angle triangle convolved with the Gaussian. The comparison uses a tolerance proportional to machine precision and the offset norm. The Gaussian and triangular kernels each have unit area, so the output area is `ampl` (including its sign).

## Inputs and units

- `offs` - one or three finite real numeric values; a three-value input specifies triangle vertices.
- `ampl` - finite real numeric scalar multiplier.
- `fwhm` - finite positive real numeric scalar Gaussian full width at half maximum.
- `x` - finite real numeric array of any dimension.

All offsets, `x`, and `fwhm` use the same coordinate units. There is no Hz-to-angular-frequency conversion; keep the caller's frequency convention consistent. For dimensionless `ampl`, `y` has reciprocal-coordinate units.

## Output

- `y` - line-shape values with the same size as `x`.
