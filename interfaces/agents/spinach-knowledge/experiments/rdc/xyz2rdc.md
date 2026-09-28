# experiments/rdc/xyz2rdc.m

- Signature: `rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,order_spec)`

## Purpose

Computes the weak heteronuclear residual dipolar coupling for two spins from their Cartesian coordinates and an order specification.

## Physical / mathematical content

For the supported `saupe` option, the routine obtains the dipolar coupling tensor `D` in rad/s using `xyz2dd` and evaluates `rdc=(2/3)*trace(S*D)/(2*pi)`, returning Hz. The Saupe matrix `S` is real, symmetric, traceless, and dimensionless.

## Parameters / inputs

- `spin_a`, `spin_b` — character strings specifying the two different isotope types (for example, `'13C'`).
- `xyz_a`, `xyz_b` — three-element Cartesian coordinate vectors for the two spins, in Angstroms.
- `order_spec` — cell array `{S,'saupe'}`, where `S` is the Saupe order matrix.

## Output

- `rdc` — weak heteronuclear residual dipolar coupling in Hz.
