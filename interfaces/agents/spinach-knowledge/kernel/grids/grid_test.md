# kernel/grids/grid_test.m

- Signature: `grid_profile=grid_test(alphas,betas,gammas,weights,ranks,sfun)`

## Purpose

Evaluates grid integration quality at each requested spherical rank using residual norms of integrated Wigner functions or spherical harmonics. If no output is requested, plots integration residual against spherical rank.

## Parameters / inputs

- `alphas` — alpha Euler angles in radians; zeros for single-angle grids.
- `betas` — beta Euler angles in radians.
- `gammas` — gamma Euler angles in radians; zeros for two-angle grids.
- `weights` — positive grid-point weights.
- `ranks` — vector of non-negative integer spherical ranks to consider.
- `sfun` — `'D_lmn'` for three-angle grids, `'Y_lm'` for two-angle grids, or `'Y_l0'` for single-angle grids.

The angle and weight inputs must be finite real column vectors of equal length.

## Output

- `grid_profile` — vector of residual norms, one per requested rank.

## Implementation

For each rank, the function sums weighted Wigner matrices over grid points. It computes the residual from the full matrix norm for `'D_lmn'`, the central-row norm for `'Y_lm'`, or the central-element norm for `'Y_l0'`, subtracting the rank-zero Kronecker delta. It reports each result and plots the profile when called without an output argument.

- Author: ilya.kuprov@weizmann.ac.il
- [Source documentation](https://spindynamics.org/wiki/index.php?title=grid_test.m)