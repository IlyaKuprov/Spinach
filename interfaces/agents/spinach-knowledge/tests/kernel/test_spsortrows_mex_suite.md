# tests/kernel/test_spsortrows_mex_suite.m

## Purpose

Regression test suite for the sparse sortrows MEX helper `spsortrows`, verifying that it returns the same row permutation as MATLAB's built-in `sortrows` on sparse real double matrices.

## Behaviour

- Announces the test target with `fprintf('TESTING: Sparse sortrows MEX helper\n')`.
- Initialises a test result via `new_test_result('kernel/spsortrows_mex_suite', 'Sparse sortrows MEX helper', 'spsortrows must return the same row permutation as Matlab sortrows.')`.
- Checks empty matrices: `sparse(0,0)` (empty square) and `sparse(0,5)` (empty tall, zero-row), comparing `spsortrows(A)` against the second output of `sortrows(A)`.
- Checks zero-column matrices with `sparse(4,0)`.
- Checks duplicate rows and lexicographic sign ordering using a 5-by-4 sparse matrix containing repeated rows, negative entries, and a zero row.
- Checks missing-value ordering with a 5-by-3 sparse matrix containing `NaN`, `Inf`, and `-Inf`, using `isequaln` for the comparison.
- Checks random sparse matrices: seeds the RNG with `rng(1729,'twister')` (saving and restoring the prior state via `onCleanup`), then runs 80 iterations. Each iteration builds `sprandn(n_rows, n_cols, density)` with `n_rows` drawn from `[1 60]`, `n_cols` from `[1 40]`, and `density = 10^(-2.5 + 2*rand)`. Every 7th iteration sets the first column to a zero sparse column; every 11th iteration sets a random element to `NaN`; every 13th sets a random element to `Inf`. The loop breaks early if any `spsortrows(A)` result differs from the `sortrows` reference permutation.
- Each check is recorded with `test_true`, with failure messages such as 'empty sparse matrices should return an empty permutation', 'zero-row sparse matrices should preserve Matlab index shape', 'zero-column sparse matrices should keep the input row order', 'duplicate sparse rows should keep Matlab stable ordering', 'NaN and Inf ordering must match Matlab sortrows', and 'random sparse real double matrices should match Matlab sortrows exactly'.

## Inputs and outputs

- Inputs: none. The function takes no arguments.
- Outputs: `result` — regression test result with explanatory messages, accumulated through `test_true` calls.

## References

- Source: [tests/kernel/test_spsortrows_mex_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_spsortrows_mex_suite.m)
