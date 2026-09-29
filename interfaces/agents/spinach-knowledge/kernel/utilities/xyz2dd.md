# kernel/utilities/xyz2dd.m

## Purpose

Converts a coordinate specification of the dipolar interaction into the dipolar interaction constant, three Euler angles, and the dipolar interaction matrix.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/xyz2dd.m>

## Behaviour

- Syntax: `[d,alp,bet,gam,M]=xyz2dd(r1,r2,isotope1,isotope2)`.
- Validates inputs via an internal consistency check (`grumble`): each coordinate vector must be numeric, real, and contain exactly three elements; `r1` and `r2` must have the same size; the two positions must differ (nonzero Euclidean distance); both isotope specifications must be character strings.
- Uses fundamental constants `hbar=1.054571628e-34` and `mu0=4*pi*1e-7`.
- Computes the inter-spin distance as the 2-norm of `r2-r1` and the unit vector `ort=(r2-r1)/distance`.
- Computes the dipolar interaction constant as `spin(isotope1)*spin(isotope2)*hbar*mu0/(4*pi*(distance*1e-10)^3)`, i.e. the distance in Angstroms is converted to meters via the `1e-10` factor.
- Derives Euler angles from the unit vector using `cart2sph`, with `bet=pi/2-bet` and `gam=0`.
- If more than four output arguments are requested, builds the 3x3 dipolar coupling matrix `M=d*[1-3*ort_i*ort_j ...]`, then symmetrises it (`M=(M+M')/2`) and removes its trace (`M=M-eye(3)*trace(M)/3`) to clean up rounding errors.
- Notes from the header: Euler angles are not uniquely defined for the orientation of axial interactions (the gamma angle can be anything); free-particle magnetogyric ratios are used, and `xyz2hfc.m` should be used instead if the system contains electrons.

## Inputs and outputs

Inputs:

- `r1`, `r2` — 3-element vectors of spin coordinates in Angstroms.
- `isotope1`, `isotope2` — isotope specification strings, e.g. `'13C'`.

Outputs:

- `d` — dipolar coupling constant, rad/s.
- `alp` — alpha Euler angle, radians.
- `bet` — beta Euler angle, radians.
- `gam` — gamma Euler angle, radians.
- `M` — dipolar interaction tensor, rad/s (computed only when requested).

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=xyz2dd.m>
- Source file: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/xyz2dd.m>
