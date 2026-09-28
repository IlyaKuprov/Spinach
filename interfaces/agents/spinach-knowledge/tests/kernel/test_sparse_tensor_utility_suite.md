# tests/kernel/test_sparse_tensor_utility_suite.m

- Signature: `result=test_sparse_tensor_utility_suite()`

## Purpose

Checks sparse, tensor, spectral-density, and SVD utilities against small explicit reference cases.

## Physical / mathematical content

The rank-dependent spectral density is checked against its Lorentzian definition using `L=2`, `Drot=1.5e6`, and `omega=2.0e5`. Blicharski invariants are checked against their antisymmetric- and traceless-symmetric-tensor definitions, and `blprod` against the corresponding polarisation identities.

## Numerical / algorithmic content

Both `kronm` and `kronm_new` are compared with an explicitly formed Kronecker product of two sparse 2-by-2 matrices, at absolute and relative tolerances of `1e-14`. The SVD truncation check uses singular values `[5 1 0.01]`: with Frobenius tolerance `0.02`, `frob_chop` must retain rank 2. For `svd_shrink`, singular values below `1e-4` are removed from `diag([4 1 1e-6])`; the returned factors must reconstruct `diag([4 1 0])` within `1e-12` absolute and relative tolerances.

## Outputs

`result` is the regression-test record with explanatory messages.

## Implementation structure

Computes each utility result from a small fixture and compares it with direct algebraic references, including the retained-singular-value reconstruction.
