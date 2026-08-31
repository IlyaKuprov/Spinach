# tests/kernel/test_spunicols_mex_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_spunicols_mex_suite.m`
- Signature: `result=test_spunicols_mex_suite()`
- Total lines: 92

## Purpose

Tests the sparse unique-column MEX helper. Syntax: result=test_spunicols_mex_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Sparse unique-column MEX helper\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/spunicols_mex_suite', 'Sparse unique-column MEX helper', 'spunicols must match Matlab unique(A.'',''rows'').'' for sparse real double matr…`.
- Lines 25-26: Check empty matrices; implemented by `A=sparse(0,0); A_ref=unique(A.','rows').'`.
- Lines 33-34: Check zero-column and all-zero matrices; implemented by `A=sparse(4,0); A_ref=unique(A.','rows').'`.
- Lines 41-42: Check duplicate columns and lexicographic signs; implemented by `A=sparse([0 0 0 0 0 0;1 1 -1 0 -1 1;0 0 3 0 3 0;-2 -2 0 0 0 -2])`.
- Lines 47-48: Check missing-value ordering and NaN duplicate retention; implemented by `A=sparse([0 0 NaN NaN 0;Inf Inf 0 0 -Inf;-Inf -Inf 1 1 1;0 0 2 2 0])`.
- Lines 53-54: Check output sparsity; implemented by `A=sprandn(30,40,0.05); A(:,20:30)=A(:,1:11)`.
- Lines 59-60: Check random sparse matrices against Matlab; implemented by `rng_state=rng`.

### Control flow inferred from the code

- Line 64: `for` loop over `n=1:120`.
- Line 69: conditional branch on `(n_cols>1)&&(mod(n,5)==0)`.
- Line 72: conditional branch on `(n_cols>2)&&(mod(n,7)==0)`.
- Line 75: conditional branch on `mod(n,11)==0`.
- Line 78: conditional branch on `mod(n,13)==0`.
- Line 82: conditional branch on `~isequaln(spunicols(A),A_ref)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/spunicols_mex_suite', 'Sparse unique-column MEX helper', 'spunicols must match Matlab unique(A.'',''rows'').'' for sparse real double matr…`.
- Lines 26: computes `A` using `A=sparse(0,0); A_ref=unique(A.','rows').'`.
- Lines 43: computes `A_ref` using `A_ref=unique(A.','rows').'`.
- Lines 55: computes `A_mex` using `A_mex=spunicols(A)`.
- Lines 60: computes `rng_state` using `rng_state=rng`.
- Lines 61: computes `rng_cleanup` using `rng_cleanup=onCleanup(@()rng(rng_state))`.
- Lines 63: computes `rand_ok` using `rand_ok=true`.
- Lines 65: computes `n_rows` using `n_rows=randi([1 80])`.
- Lines 66: computes `n_cols` using `n_cols=randi([1 70])`.
- Lines 67: computes `density` using `density=10^(-2.7+2.2*rand)`.
- Lines 70: computes `A(:,n_cols)` using `A(:,n_cols)=A(:,1)`.
- Lines 73: computes `A(:,2)` using `A(:,2)=sparse(n_rows,1)`.
- Lines 76: computes `A(randi(n_rows),randi(n_cols))` using `A(randi(n_rows),randi(n_cols))=NaN`.

## Outputs

- result -regression test result with explanatory messages
- The test compares spunicols() against Matlab unique(A.','rows').' on
- empty, zero-column, duplicate-column, NaN, Inf, signed-value, and random
- sparse real double matrices.

## Implementation structure

- Tests the sparse unique-column MEX helper. Syntax:
- result=test_spunicols_mex_suite()
- result -regression test result with explanatory messages
- The test compares spunicols() against Matlab unique(A.','rows').' on
- empty, zero-column, duplicate-column, NaN, Inf, signed-value, and random
- sparse real double matrices.
- Announce the test target
- State the utility target of the test
- Check empty matrices
- Check zero-column and all-zero matrices
- Check duplicate columns and lexicographic signs
- Check missing-value ordering and NaN duplicate retention

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_true()`, `isequal()`, `spunicols()`, `isequaln()`, `sprandn()`, `issparse()`, `onCleanup()`, `rng()`, `randi()`, `clear()`.
