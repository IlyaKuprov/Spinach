# kernel/conventions/transforms/anas2mat.m

- Signature: `M = anas2mat(iso,an,as,alp,bet,gam)`
- Source: [`kernel/conventions/transforms/anas2mat.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/anas2mat.m)
- Existing Wiki page: [`anas2mat.m`](https://spindynamics.org/wiki/index.php?title=anas2mat.m)

## Contract

This utility turns the isotropic part, anisotropy, asymmetry, and Euler orientation of a real 3x3 interaction tensor into its Cartesian matrix. The input components are scalar values in the same units as the represented interaction; the asymmetry is dimensionless. Euler angles are in radians, and the orientation convention is the one implemented by `euler2dcm`.

The principal values are reconstructed using `ra = 2*an/3`, then `zz = iso+ra`, `yy = iso-ra*(1-as)/2`, and `xx = iso-ra*(1+as)/2`. Thus `iso = (xx+yy+zz)/3`, `an = zz-(xx+yy)/2`, and `as = (yy-xx)/(zz-iso)`, in the Haeberlen-Mehring convention used by the source. It then computes `R = euler2dcm(alp,bet,gam)` and returns `M = R*diag([xx yy zz])*R'`.

## Inputs and output

- `iso`: isotropic tensor component, the mean of the three principal values.
- `an`: anisotropy, defined from the principal values as above.
- `as`: dimensionless asymmetry parameter.
- `alp`, `bet`, `gam`: alpha, beta, and gamma Euler angles in radians.
- `M`: real 3x3 tensor matrix in the Cartesian frame produced by the stated rotation.

All six inputs must be real numeric scalars. The function performs this shape/type check but does not impose a particular unit system or additional physical range for the asymmetry.

## Source-supported example

At zero Euler angles, the rotation is the identity, so `M = diag([xx yy zz])` with `xx`, `yy`, and `zz` given by the reconstruction equations above. This illustrates the principal-axis input without choosing an interaction unit or a special asymmetry value.
