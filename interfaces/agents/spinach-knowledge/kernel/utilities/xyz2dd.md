# kernel/utilities/xyz2dd.m

- Signature: `[d,alp,bet,gam,M]=xyz2dd(r1,r2,isotope1,isotope2)`

Converts two spin coordinates and their isotope specifications into a dipolar coupling constant, three Euler angles, and, if requested, a dipolar interaction tensor.

## Inputs

- `r1`, `r2`: Three-element real coordinate vectors of the same dimensions, in ångströms. The coordinates must differ.
- `isotope1`, `isotope2`: Isotope specification character strings, such as `'13C'`.

## Outputs

- `d`: Dipolar coupling constant, in rad/s.
- `alp`, `bet`, `gam`: Euler angles, in radians. For this axial interaction, the angles are not unique; the function sets `gam=0`.
- `M`: Dipolar interaction tensor, in rad/s, computed when requested as a fifth output.

## Calculation

The function uses the separation `distance=norm(r2-r1,2)` and unit direction `ort=(r2-r1)/distance`. It computes `d` from the product of the isotope magnetogyric ratios returned by `spin`, `hbar`, `mu0`, and the inverse cube of the separation converted from ångströms to metres. The direction determines `alp` and `bet`; `gam` is set to zero. When requested, the tensor is `M=d*(eye(3)-3*ort(:)*ort(:)')`, then symmetrized and made traceless to remove rounding errors.

Free-particle magnetogyric ratios are used. For systems containing electrons, use `xyz2hfc.m` instead.

[Source page](https://spindynamics.org/wiki/index.php?title=xyz2dd.m)