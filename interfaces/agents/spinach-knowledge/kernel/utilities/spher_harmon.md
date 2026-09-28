# kernel/utilities/spher_harmon.m

- Signature: `Y=spher_harmon(l,m,theta,phi)`

## Purpose

Evaluate spherical harmonics at the specified angles.

## Physical / mathematical content

- Uses Schmidt-normalized associated Legendre functions and the azimuthal factor `exp(1i*m*phi)`.

## Numerical / algorithmic content

- Computes `S=legendre(l,cos(theta),'sch')` and selects the `abs(m)+1` component, reshaping it to the size of `theta`.
- Forms `Y=sqrt((2*l+1)/(4*pi))*S.*exp(1i*m*phi)`, dividing by `sqrt(2)` when `m` is nonzero. Negates `Y` when `m` is positive and odd.

## Parameters / inputs

- `l` — L quantum number; a numeric, real, scalar, non-negative integer.
- `m` — M quantum number; a numeric, real, scalar integer in `[-l,l]`.
- `theta` — numeric, real array of theta angles in radians.
- `phi` — numeric, real array of phi angles in radians.

## Outputs

- `Y` — array of spherical harmonics evaluated at the specified angles.

## Implementation structure

- Checks input constraints with `grumble(l,m,theta,phi)`, computes the selected Schmidt-normalized Legendre component, then forms the spherical harmonics.

Source: [spher_harmon.m](https://spindynamics.org/wiki/index.php?title=spher_harmon.m).

ilya.kuprov@weizmann.ac.il

<https://spindynamics.org/wiki/index.php?title=spher_harmon.m>