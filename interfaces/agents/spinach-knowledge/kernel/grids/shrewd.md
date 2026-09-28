# kernel/grids/shrewd.m

- Signature: `weights=shrewd(alphas,betas,gammas,max_rank,max_error)`

## Purpose

Computes SHREWD weights for a two- or three-angle spherical grid. For the algorithm, see the paper by Eden and Levitt: http://dx.doi.org/10.1006/jmre.1998.1427. Function page: https://spindynamics.org/wiki/index.php?title=shrewd.m.

## Inputs

- `alphas`, `betas`, `gammas`: matching column vectors of finite, real Euler angles in radians, using the active ZYZ convention. Set `alphas` to all zeros for a two-angle grid.
- `max_rank`: finite positive integer giving the maximum spherical rank considered when minimizing residuals.
- `max_error`: finite non-negative real scalar giving the maximum residual absolute error per spherical function.

## Output

- `weights`: one grid weight for each supplied `[alpha beta gamma]` point.

## Implementation

The function validates its inputs, then selects the two-angle branch when every value in `alphas` is zero. That branch assembles a complex spherical-harmonic matrix from Wigner-function entries through `max_rank`. Otherwise, it assembles a Wigner-function matrix over rank and both magnetic indices. In either branch, the right-hand-side vector contains `max_error` in every position except the first, which contains `1-max_error`.

The function solves the matrix system for the weights, takes their real parts, and normalizes them to sum to one. It raises an error if any resulting weight is zero, advising an increase in maximum rank, or negative, advising a reduction in the accuracy threshold.