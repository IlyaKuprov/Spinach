# tests/kernel/test_spunicols_mex_suite.m

- Signature: `result=test_spunicols_mex_suite()`

## Purpose

Tests the sparse unique-column MEX helper `spunicols` against MATLAB's `unique(A.','rows').'` for sparse real double matrices.

## Numerical / algorithmic content

The tests cover empty and zero-column matrices, all-zero and duplicate columns, signed values, NaN and Inf ordering, output sparsity, and randomized sparse matrices. Comparisons use `isequal` or `isequaln` where NaNs may occur.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announces the test target and initializes the regression result.
- Checks empty, zero-column, all-zero, duplicate-column, and nonfinite-value cases against the MATLAB reference.
- Verifies that `spunicols` returns a sparse matrix.
- Runs 120 seeded random sparse-matrix comparisons, including injected duplicate or zero columns and occasional NaN or Inf values; restores the prior random-number generator state afterward.