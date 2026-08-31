# experiments/rdc/rdc_fit.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/rdc/rdc_fit.m`
- Signature: `S=rdc_fit(isotopes,xyz,rdc)`
- Total lines: 88

## Purpose

Linear least squares fitter for residual dipolar couplings. Iso- tope pairs are arbitrary heteronuclear; multiple isotope pairs may be supplied at the same time. Syntax: S=rdc_fit(isotopes,xyz,rdc)

## Physical / mathematical content

- Residual-dipolar-coupling experiment and analysis routines. These files use partial ordering, Saupe tensors, and molecular-frame geometry to connect internuclear vectors with observed couplings.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(isotopes,xyz,rdc)`.
- Lines 32-33: Build dipolar tensor array; implemented by `D=cell(1,size(isotopes,1))`.
- Lines 36-37: Get the dipole-dipole coupling tensor; implemented by `[~,~,~,~,D{n}]=xyz2dd(xyz{n,1},xyz{n,2},isotopes{n,1},isotopes{n,2})`.
- Lines 39-40: Extract second rank component; implemented by `[~,~,D{n}]=mat2sphten(D{n})`.
- Lines 44-45: Run pseudoinverse least squares; implemented by `S=(3/2)*(cell2mat(D)'\(2*pi*rdc(:)))`.
- Lines 47-48: Convert back into a matrix; implemented by `S=real(sphten2mat([],[],S))`.

### Control flow inferred from the code

- Line 34: `for` loop over `n=1:size(isotopes,1)`.

### Key state/data transformations

- Lines 33: computes `D` using `D=cell(1,size(isotopes,1))`.
- Lines 37: computes `[~,~,~,~,D{n}]` using `[~,~,~,~,D{n}]=xyz2dd(xyz{n,1},xyz{n,2},isotopes{n,1},isotopes{n,2})`.
- Lines 40: computes `[~,~,D{n}]` using `[~,~,D{n}]=mat2sphten(D{n})`.
- Lines 45: computes `S` using `S=(3/2)*(cell2mat(D)'\(2*pi*rdc(:)))`.

### Local helper functions

- Line 53: `grumble()` — `function grumble(isotopes,xyz,rdc)`.
  - Representative operation: `if (~iscell(isotopes))||(size(isotopes,2)~=2)`.
  - Representative operation: `error('isotopes must be an N x 2 cell array of character strings.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xyz2dd()`, `mat2sphten()`, `cell2mat()`, `rdc()`, `sphten2mat()`, `iscell()`, `ischar()`, `strcmp()`, `iscolumn()`.
