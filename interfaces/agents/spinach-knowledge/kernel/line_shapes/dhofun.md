# kernel/line_shapes/dhofun.m

- MATLAB source: [kernel/line_shapes/dhofun.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/line_shapes/dhofun.m)
- Existing Wiki: [dhofun.m](https://spindynamics.org/wiki/index.php?title=dhofun.m)
- Signature: `y=dhofun(x,nat_freq,fwhm)`

## Meaning and equation

This is the normalised damped-harmonic-oscillator response in magnetic-resonance notation. For `x>0`, put `r=x/nat_freq` and `d=fwhm/nat_freq`; the source evaluates `y=(2*d/(pi*nat_freq))*r^2/((r^2-1)^2+(d*r)^2)`. It sets `y=0` for non-positive `x`. The response integrates to one over positive arguments, peaks at `nat_freq`, and the source identifies `fwhm` as the full width at half maximum at any damping; in the weak-damping limit it tends to a Lorentzian of that width.

## Inputs and units

- `x` - real numeric array of any dimension. The guard checks that it is numeric and real; it does not explicitly reject non-finite values.
- `nat_freq` - finite positive real numeric scalar.
- `fwhm` - finite positive real numeric scalar.

`x`, `nat_freq`, and `fwhm` must use the same frequency coordinate. The function does not convert between Hz and angular frequency, so the caller's consistent convention determines which is used. The response has reciprocal-frequency units.

## Output

- `y` - same size and type as `x`, initialised with zeros like `x`; positive-argument entries are replaced by the response.
