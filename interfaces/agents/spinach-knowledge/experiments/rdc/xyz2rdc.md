# experiments/rdc/xyz2rdc.m

- Signature: `rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,order_spec)`

## Purpose

Converts Cartesian coordinates of a pair of nuclei and an order matrix into residual dipolar coupling; only hetero- nuclear spin pairs are supported. Syntax: rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,chi)

## Physical / mathematical content

- Residual-dipolar-coupling experiment and analysis routines. These files use partial ordering, Saupe tensors, and molecular-frame geometry to connect internuclear vectors with observed couplings.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

## Parameters / inputs

- spin_a, spin_b -character strings indicating spin
- type, for example '13C'
- xyz_a, xyz_b -three-element vectors specifying
- Cartesian coordinates of the two
- spins in Angstroms
- order_spec -{S,'saupe'} uses Saupe order mat-
- rix, S is a traceless symmetric
- 3x3 matrix, dimensionless

## Outputs

- rdc -residual dipolar coupling in the
- heteronuclear case, Hz

## Implementation structure

- Converts Cartesian coordinates of a pair of nuclei and an
- order matrix into residual dipolar coupling; only hetero-
- nuclear spin pairs are supported. Syntax:
- rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,chi)
- spin_a, spin_b -character strings indicating spin
- type, for example '13C'
- xyz_a, xyz_b -three-element vectors specifying
- Cartesian coordinates of the two
- spins in Angstroms
- order_spec -{S,'saupe'} uses Saupe order mat-
- rix, S is a traceless symmetric
- 3x3 matrix, dimensionless
