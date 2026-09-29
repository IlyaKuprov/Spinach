# experiments/rdc/rdc_fit.m

Source: [experiments/rdc/rdc_fit.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/rdc/rdc_fit.m)
Spinach Wiki: [rdc_fit.m](https://spindynamics.org/wiki/index.php?title=rdc_fit.m)

- Signature: `S=rdc_fit(isotopes,xyz,rdc)`

## Purpose

Fits one Saupe order matrix to a set of measured residual dipolar couplings by linear least squares. Each row may use a different heteronuclear isotope pair; this is a tensor fit to supplied geometry and RDC values, not a pulse-sequence calculation.

## Inputs and coordinate convention

- `isotopes` is an N-by-2 cell array of isotope-name strings (for example, `'13C'` labels).
- `xyz` is an N-by-2 cell array; each entry is a real three-element Cartesian coordinate vector in Angstroms. The coordinates for all pairs must use a common Cartesian frame, which is the frame of the fitted order matrix.
- `rdc` is an N-by-1 real column of measured couplings in Hz.

The three inputs must have matching row counts. Isotope pairs are heteronuclear; the pair may vary between rows.

## Fit

For each row, the routine obtains the dipole-dipole coupling tensor from `xyz2dd` and takes its rank-2 spherical-tensor components with `mat2sphten`. It assembles those components into a linear design matrix and solves with MATLAB backslash using the scaled observations `2*pi*rdc`; the fitted rank-2 coefficients are converted to a matrix with `sphten2mat`. In code terms, the solve is `(3/2)*(cell2mat(D)' \ (2*pi*rdc(:)))` before that conversion.

The output `S` is a real, dimensionless 3-by-3 Saupe order matrix represented by the rank-2 tensor basis. The implementation performs a linear fit without an explicit regularizer or nonlinear refinement.
