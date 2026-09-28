# kernel/conventions/transforms/xyz2sph.m

- Signature: `[r,theta,phi] = xyz2sph(x,y,z)`

Converts Cartesian coordinates to spherical coordinates using the ISO convention.

## Inputs

- `x`, `y`, `z`: arrays of Cartesian X, Y, and Z coordinates. All must be numeric, real, and the same size.

## Outputs and conversion

- `r`: radius, `sqrt(x.^2+y.^2+z.^2)`; documented range `0 <= r < Inf`.
- `theta`: inclination, `acos(z./r)`; documented range `0 <= theta <= pi`.
- `phi`: azimuth, `mod(atan2(y,x),2*pi)`; range `0 <= phi < 2*pi`.

The formulas operate elementwise. At the origin, `theta` evaluates to `NaN` because `r` is zero.

[Source](https://spindynamics.org/wiki/index.php?title=xyz2sph.m)