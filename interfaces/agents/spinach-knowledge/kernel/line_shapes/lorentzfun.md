# kernel/line_shapes/lorentzfun.m

- Signature: `[real_part,imag_part]=lorentzfun(offs,ampl,fwhm,x,phi)`

## Meaning

Returns the absorption and dispersion components of a Lorentzian line with a phase rotation. Let `gamma=fwhm/2` and `u=(x-offs)/gamma`. The unrotated Lorentzian is `L=ampl/(2*pi*gamma)/(1+u^2)`; the returned arrays are `real_part=L*cos(phi)-u*L*sin(phi)` and `imag_part=L*sin(phi)+u*L*cos(phi)`. At zero phase the absorption component integrates to `ampl/2`, as stated in the source's one-sided FID Fourier-transform convention; the source contrasts this with `lorentzcon()`, whose integral is `ampl`.

## Inputs and outputs

- `offs`: real numeric scalar peak offset; `ampl`: real numeric scalar amplitude multiplier.
- `fwhm`: positive real numeric scalar full width at half maximum.
- `x`: real numeric array of any dimension. Both outputs have the same size as `x`.
- `phi`: real numeric scalar phase in radians.

The formula performs no conversion between Hz and angular frequency. Use the same coordinate units for `offs`, `fwhm`, and `x`; the source does not specify which frequency convention is intended. The checks enforce the real/numeric/scalar or array conditions above and positive width.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/line_shapes/lorentzfun.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=lorentzfun.m)
