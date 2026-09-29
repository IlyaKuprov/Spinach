# experiments/rdc/xyz2rdc.m

Source: [experiments/rdc/xyz2rdc.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/rdc/xyz2rdc.m)
Spinach Wiki: [xyz2rdc.m](https://spindynamics.org/wiki/index.php?title=xyz2rdc.m)

- Signature: `rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,order_spec)`

## Purpose

Calculates one weak heteronuclear residual dipolar coupling from a nuclear-spin pair's Cartesian coordinates and a Saupe order matrix. This is a geometry-to-coupling utility, not an RDC fitting routine.

## Inputs and coordinate convention

- `spin_a` and `spin_b` are isotope-name strings and must identify different isotopes (for example, one may be `'13C'`).
- `xyz_a` and `xyz_b` are real three-element coordinate vectors in Angstroms.
- `order_spec` is a cell specification `{S,'saupe'}`. The consistency check accepts a 2- or 4-element cell and the calculation uses its first two entries. The documented Saupe matrix `S` is dimensionless, symmetric, traceless, and 3-by-3. Coordinates and `S` must be expressed in the same Cartesian frame. The implementation accepts a real 3-by-3 matrix but does not test its symmetry or trace.

## Calculation and output

The routine obtains the dipolar coupling tensor `D` from `xyz2dd` (rad/s), then evaluates `rdc=(2/3)*trace(S*D)/(2*pi)`. The result is the heteronuclear residual dipolar coupling in Hz. The only supported order-specification type is `'saupe'`.
