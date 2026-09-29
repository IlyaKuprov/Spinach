# kernel/conventions/transforms/xyz2sph.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/xyz2sph.m) · [Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=xyz2sph.m)

## Purpose and equations

Converts Cartesian coordinates to spherical coordinates in the ISO convention, element by element:

- `r=sqrt(x.^2+y.^2+z.^2)`
- `theta=acos(z./r)` (inclination from the positive Z axis; nominal range `0<=theta<=pi`)
- `phi=mod(atan2(y,x),2*pi)` (azimuth from the positive X axis toward positive Y; range `0<=phi<2*pi`)

`x`, `y`, and `z` must be numeric, real arrays with identical sizes. Outputs have the same array shape. The radius has the same units as the input coordinates; the angles are in radians, with no unit conversion applied. The implementation does not check finiteness. At the origin, `r=0` and `theta` evaluates as `NaN` from `0/0`; azimuth is not geometrically defined there even though the formula is evaluated.
