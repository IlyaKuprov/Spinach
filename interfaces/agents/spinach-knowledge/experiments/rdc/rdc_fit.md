# experiments/rdc/rdc_fit.m

- Signature: `S=rdc_fit(isotopes,xyz,rdc)`

## Purpose

Linear least squares fitter for residual dipolar couplings. Iso- tope pairs are arbitrary heteronuclear; multiple isotope pairs may be supplied at the same time. Syntax: S=rdc_fit(isotopes,xyz,rdc)

## Physical / mathematical content

- Residual-dipolar-coupling experiment and analysis routines. These files use partial ordering, Saupe tensors, and molecular-frame geometry to connect internuclear vectors with observed couplings.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

## Parameters / inputs

- isotopes -N x 2 cell array of strings with Spinach
- isotope specifications, e.g. '13C'
- xyz -N x 2 cell array of 3-element vectors
- with Cartesian coordinates in Angstrom
- rdc -N x 1 vector with residual dipolar coup-
- lings in Hz

## Outputs

- S -Saupe order matrix, a symmetric real
- dimensionless 3x3 matrix

## Implementation structure

- Linear least squares fitter for residual dipolar couplings. Iso-
- tope pairs are arbitrary heteronuclear; multiple isotope pairs
- may be supplied at the same time. Syntax:
- S=rdc_fit(isotopes,xyz,rdc)
- isotopes -N x 2 cell array of strings with Spinach
- isotope specifications, e.g. '13C'
- xyz -N x 2 cell array of 3-element vectors
- with Cartesian coordinates in Angstrom
- rdc -N x 1 vector with residual dipolar coup-
- lings in Hz
- S -Saupe order matrix, a symmetric real
- dimensionless 3x3 matrix
