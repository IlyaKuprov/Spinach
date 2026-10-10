# kernel/line_shapes/lorentzcon.m

- MATLAB source: [kernel/line_shapes/lorentzcon.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/line_shapes/lorentzcon.m)
- Existing Wiki: [lorentzcon.m](https://spindynamics.org/wiki/index.php?title=lorentzcon.m)
- Signature: `y=lorentzcon(offs,ampl,fwhm,x)`

## Meaning and equation

The routine returns a Lorentzian, or its convolution with a normalised boxcar or triangular distribution. Let `gamma=fwhm/2`; the unit-area kernel centred at zero is `L(u)=gamma/(pi*(u^2+gamma^2))`.

- One offset `q`: `y(x)=ampl*L(x-q)`.
- Two sorted offsets `a<b`: they define a uniform density `B(t)=1/(b-a)` on `[a,b]`; `y(x)=ampl*integral(B(t)*L(x-t),t=a..b)`. This is the source's Lorentzian-boxcar convolution.
- Three distinct sorted offsets `a<b<c`: they define the unit-area triangular density `T(t)=2*(t-a)/((b-a)*(c-a))` on `[a,b]`, `T(t)=2*(c-t)/((c-b)*(c-a))` on `(b,c]`, and zero elsewhere; `y(x)=ampl*integral(T(t)*L(x-t),t=a..c)`.

For three offsets, repeated or numerically coalescent vertices select the limiting Lorentzian or right-angle-triangle convolution. The implementation compares spacings with a tolerance based on machine precision and the offset norm. Since each kernel is unit area, the line-shape area is `ampl`.

## Inputs and units

- `offs` - one, two, or three finite real numeric values.
- `ampl` - finite real numeric scalar multiplier.
- `fwhm` - finite positive real numeric scalar full width at half maximum.
- `x` - finite real numeric array of any dimension.

The offsets, `x`, and `fwhm` must share the same coordinate units. The function does not convert between Hz and angular frequency; choose and use one convention consistently. For dimensionless `ampl`, `y` has reciprocal-coordinate units. The source converts the offsets, amplitude, and width to double precision and converts integer `x` to double before evaluation.

## Output

- `y` - values with the same size as `x` (integer `x` is promoted to double by the source).
