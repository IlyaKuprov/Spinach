# tests/kernel/test_spsortrows_mex_suite.m

- Signature: `result=test_spsortrows_mex_suite()`

## Purpose

Regression-tests the sparse `spsortrows` MEX helper by comparing its returned row permutation with MATLAB `sortrows` for sparse real double matrices.

## Numerical / algorithmic content

The suite checks empty `0`-by-`0` and `0`-by-`5` matrices, a `4`-by-`0` matrix, duplicate rows with positive and negative values, and rows containing `NaN`, `Inf`, and `-Inf`. It also runs 80 reproducible random sparse-matrix comparisons. Random matrices have 1–60 rows, 1–40 columns, and density `10^(-2.5+2*rand)`; selected iterations add an all-zero column, `NaN`, or `Inf`. The test saves and restores the MATLAB random-number-generator state and uses seed `1729` with the `twister` generator.

## Outputs

- `result` — regression test result with explanatory messages. Each check compares `spsortrows(A)` with the row indices returned by `sortrows(A)`; the nonfinite-value check uses `isequaln`, and the other checks use `isequal`. The random-set check requires every comparison to match.