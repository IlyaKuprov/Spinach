# tests/kernel/test_spsortrows_mex_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_spsortrows_mex_suite.m`
- Signature: `result=test_spsortrows_mex_suite()`
- Total lines: 84

## Purpose

Tests the sparse sortrows MEX helper. Syntax: result=test_spsortrows_mex_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Sparse sortrows MEX helper\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/spsortrows_mex_suite', 'Sparse sortrows MEX helper', 'spsortrows must return the same row permutation as Matlab sortrows.')`.
- Lines 25-26: Check empty matrices; implemented by `A=sparse(0,0)`.
- Lines 35-36: Check zero-column matrices; implemented by `A=sparse(4,0)`.
- Lines 41-42: Check duplicate rows and lexicographic signs; implemented by `A=sparse([0 1 0 -2;0 1 0 -2;0 -1 3 0;0 -1 2 9;0 0 0 0])`.
- Lines 47-48: Check missing-value ordering; implemented by `A=sparse([0 NaN 0;0 Inf 0;0 -Inf 1;0 0 2;0 NaN -1])`.
- Lines 53-54: Check random sparse matrices against Matlab; implemented by `rng_state=rng`.

### Control flow inferred from the code

- Line 58: `for` loop over `n=1:80`.
- Line 63: conditional branch on `mod(n,7)==0`.
- Line 66: conditional branch on `mod(n,11)==0`.
- Line 69: conditional branch on `mod(n,13)==0`.
- Line 73: conditional branch on `~isequal(spsortrows(A),idx_ref)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/spsortrows_mex_suite', 'Sparse sortrows MEX helper', 'spsortrows must return the same row permutation as Matlab sortrows.')`.
- Lines 26: computes `A` using `A=sparse(0,0)`.
- Lines 27: computes `[~,idx_ref]` using `[~,idx_ref]=sortrows(A)`.
- Lines 54: computes `rng_state` using `rng_state=rng`.
- Lines 55: computes `rng_cleanup` using `rng_cleanup=onCleanup(@()rng(rng_state))`.
- Lines 57: computes `rand_ok` using `rand_ok=true`.
- Lines 59: computes `n_rows` using `n_rows=randi([1 60])`.
- Lines 60: computes `n_cols` using `n_cols=randi([1 40])`.
- Lines 61: computes `density` using `density=10^(-2.5+2*rand)`.
- Lines 64: computes `A(:,1)` using `A(:,1)=sparse(n_rows,1)`.
- Lines 67: computes `A(randi(n_rows),randi(n_cols))` using `A(randi(n_rows),randi(n_cols))=NaN`.

## Outputs

- result -regression test result with explanatory messages
- The test compares spsortrows() against Matlab sortrows() on empty,
- zero-column, duplicate-row, NaN, Inf, signed-value, and random sparse
- real double matrices.

## Implementation structure

- Tests the sparse sortrows MEX helper. Syntax:
- result=test_spsortrows_mex_suite()
- result -regression test result with explanatory messages
- The test compares spsortrows() against Matlab sortrows() on empty,
- zero-column, duplicate-row, NaN, Inf, signed-value, and random sparse
- real double matrices.
- Announce the test target
- State the utility target of the test
- Check empty matrices
- Check zero-column matrices
- Check duplicate rows and lexicographic signs
- Check missing-value ordering

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `sortrows()`, `test_true()`, `isequal()`, `spsortrows()`, `isequaln()`, `onCleanup()`, `rng()`, `randi()`, `sprandn()`, `clear()`.
