# tests/kernel/test_spunicols_mex_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_spunicols_mex_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_spunicols_mex_suite.m)

## Purpose

Regression test suite for the sparse unique-column MEX helper `spunicols`. The suite verifies that `spunicols` reproduces the behaviour of the Matlab expression `unique(A.','rows').'` on sparse real double matrices, covering empty, zero-column, duplicate-column, NaN, Inf, signed-value, and random cases.

## Behaviour

- Announces the test target with `fprintf('TESTING: Sparse unique-column MEX helper\n')`.
- Initialises a test result object via `new_test_result('kernel/spunicols_mex_suite', 'Sparse unique-column MEX helper', 'spunicols must match Matlab unique(A.'',''rows'').'' for sparse real double matrices.')`.
- Empty matrices: `sparse(0,0)` and `sparse(0,5)` are compared against `unique(A.','rows').'` using `isequal`; the zero-row case expects all empty columns to collapse to one column.
- Zero-column and all-zero matrices: `sparse(4,0)` should remain a zero-column matrix, and `sparse(4,3)` (all-zero columns) should collapse to one sparse column.
- Duplicate columns and lexicographic signs: a 4-by-6 sparse matrix built from `[0 0 0 0 0 0;1 1 -1 0 -1 1;0 0 3 0 3 0;-2 -2 0 0 0 -2]` is checked with `isequal` against the Matlab reference.
- Nonfinite values: a 4-by-5 sparse matrix built from `[0 0 NaN NaN 0;Inf Inf 0 0 -Inf;-Inf -Inf 1 1 1;0 0 2 2 0]` is compared with `isequaln`, so NaN and Inf ordering must match the Matlab double-transpose reference.
- Output sparsity: for `sprandn(30,40,0.05)` with columns 20:30 overwritten by columns 1:11, the output of `spunicols` must satisfy `issparse`.
- Random sparse set: the previous `rng` state is saved and restored through `onCleanup`; the generator is seeded with `rng(1729,'twister')`. A loop of 120 iterations builds `sprandn(n_rows,n_cols,density)` matrices with `n_rows` from `randi([1 80])`, `n_cols` from `randi([1 70])`, and `density = 10^(-2.7+2.2*rand)`. Conditional perturbations: when `n_cols>1` and `mod(n,5)==0`, the last column is set equal to the first; when `n_cols>2` and `mod(n,7)==0`, column 2 is set to `sparse(n_rows,1)`; when `mod(n,11)==0`, a random element is set to NaN; when `mod(n,13)==0`, a random element is set to Inf. Each matrix is compared against `unique(A.','rows').'` using `isequaln`; the loop breaks on the first mismatch and the cleanup object is cleared afterwards.
- Each check is recorded with `test_true(result, name, condition, message)`, accumulating pass/fail outcomes and explanatory messages in the returned result.

## Inputs and outputs

```matlab
result = test_spunicols_mex_suite()
```

- **Outputs**
  - `result` — regression test result with explanatory messages.
- **Inputs**
  - None.

## References

- `spunicols` — sparse unique-column MEX helper under test.
- `new_test_result`, `test_true` — test harness utilities used to build and record the regression result.
- Matlab functions used as the reference behaviour: `unique(A.','rows').'`, `isequal`, `isequaln`, `issparse`, `sprandn`, `rng`, `randi`, `rand`, `onCleanup`.
