# kernel/utilities/xyz2dd.m

- Signature: `[d,alp,bet,gam,M]=xyz2dd(r1,r2,isotope1,isotope2)`

## Purpose

Converts coordinate specification of the dipolar interaction into the dipolar interaction constant, three Euler angles, and the dipolar interaction matrix. Syntax: [d,alp,bet,gam,m]=xyz2dd(r1,r2,isotope1,isotope2)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

## Parameters / inputs

- r1,2 -3-element vectors of spin coordinates
- in Angstroms
- isotope1,2 -isotope specification strings, e.g.
- '13C'.

## Outputs

- d -dipolar coupling constant, rad/s
- alp -alpha Euler angle, radians
- bet -beta Euler angle, radians
- gam -gamma Euler angle, radians
- M -dipolar interaction tensor, rad/s
- N.B. Euler angles are not uniquely defined for the orientati-
- on of axial interactions (gamma angle can be anything).
- N.B. free-particle magnetogyric ratios are used, use xyz2hfc.m
- if you have electrons in the system.

## Implementation structure

- Converts coordinate specification of the dipolar interaction
- into the dipolar interaction constant, three Euler angles,
- and the dipolar interaction matrix. Syntax:
- [d,alp,bet,gam,m]=xyz2dd(r1,r2,isotope1,isotope2)
- r1,2 -3-element vectors of spin coordinates
- in Angstroms
- isotope1,2 -isotope specification strings, e.g.
- '13C'.
- d -dipolar coupling constant, rad/s
- alp -alpha Euler angle, radians
- bet -beta Euler angle, radians
- gam -gamma Euler angle, radians
