# experiments/rdc/rdc_fit.m

- Signature: `S=rdc_fit(isotopes,xyz,rdc)`

## Purpose

Fits a Saupe order matrix to residual dipolar coupling (RDC) measurements by linear least squares. Heteronuclear isotope pairs may differ between observations.

## Numerical / algorithmic content

For each pair, the routine obtains the dipolar coupling tensor from `xyz2dd` and extracts its rank-2 spherical-tensor components. It solves the linear least-squares system with the `2*pi*rdc` data vector, scales the solution by `3/2`, and converts the fitted components to a real `3x3` matrix.

## Parameters / inputs

- `isotopes` — `N x 2` cell array of Spinach isotope-specification strings (for example, `'13C'`); the two isotopes in each pair must differ.
- `xyz` — `N x 2` cell array containing the two Cartesian coordinate vectors for each pair, in Angstroms.
- `rdc` — `N x 1` vector of residual dipolar couplings in Hz.

## Output

- `S` — real, symmetric, dimensionless `3x3` Saupe order matrix.
