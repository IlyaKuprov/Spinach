# kernel/conventions/transforms/sphten2mat.m

- Signature: `M=sphten2mat(rank0,rank1,rank2)`

## Purpose

Converts nine irreducible spherical-tensor components of an interaction tensor into a `3x3` Cartesian matrix. The conventions match Equation (18) of the paper by Len Mueller (http://dx.doi.org/10.1002/cmr.a.20224). Supply coefficients of the corresponding irreducible spherical-tensor operators returned by `irr_sph_ten.m` in this order:

- Rank 0: `(0,0)`
- Rank 1: `(1,1)`, `(1,0)`, `(1,-1)`
- Rank 2: `(2,2)`, `(2,1)`, `(2,0)`, `(2,-1)`, `(2,-2)`

## Parameters / inputs

- `rank0`: A single number giving the coefficient of `T(0,0)`; an empty numeric input omits this contribution.
- `rank1`: A three-element numeric vector giving the coefficients of `T(1,1)`, `T(1,0)`, and `T(1,-1)` in that order; an empty numeric input omits these contributions.
- `rank2`: A five-element numeric vector giving the coefficients of `T(2,2)`, `T(2,1)`, `T(2,0)`, `T(2,-1)`, and `T(2,-2)` in that order; an empty numeric input omits these contributions.

## Output

- `M`: `3x3` interaction tensor.

## Implementation

The function checks that all three inputs are numeric and that each is either empty or has the required number of elements (`1`, `3`, and `5`, respectively). It initializes `M` to a `3x3` zero matrix, adds `rank0*eye(3)` when `rank0` is nonempty, and adds the rank-1 and rank-2 contributions using the matrices specified in the source.

Contact: ilya.kuprov@weizmann.ac.il

Reference page: https://spindynamics.org/wiki/index.php?title=sphten2mat.m