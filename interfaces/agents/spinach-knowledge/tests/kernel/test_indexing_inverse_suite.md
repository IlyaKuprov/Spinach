# tests/kernel/test_indexing_inverse_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_indexing_inverse_suite.m`
- Signature: `result=test_indexing_inverse_suite()`
- Total lines: 102

## Purpose

Tests indexing helper inverses. Syntax: result=test_indexing_inverse_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Indexing inverse helpers\n')`.
- Lines 19-22: State the indexing target of the test; implemented by `result=new_test_result('kernel/indexing_inverse_suite', 'Indexing inverse helpers', 'indexing helpers must be exact inverses on their finite integer domains.')`.
- Lines 24-25: Check the documented base-one serpentine matrix; implemented by `S_ref=[1 3 6 10;2 5 9 13;4 8 12 15;7 11 14 16]`.
- Lines 30-31: Check the documented base-zero serpentine matrix; implemented by `S=serpentine(4,0)`.
- Lines 35-36: Check k,q to linear and back in base-one indexing; implemented by `[K,Q]=ndgrid(1:4,1:4)`.
- Lines 46-47: Check k,q to linear and back in base-zero indexing; implemented by `[K,Q]=ndgrid(0:3,0:3)`.
- Lines 57-58: Build a complete low-rank L,M domain in documented order; implemented by `L=[]; M=[]`.
- Lines 66-67: Check L,M linear indexing and inverse conversion; implemented by `I=lm2lin(L,M)`.
- Lines 76-77: Build a complete low-rank L,M,N Wigner-function domain in documented order; implemented by `L=[]; M=[]; N=[]`.
- Lines 88-89: Check L,M,N linear indexing and inverse conversion; implemented by `I=lmn2lin(L,M,N)`.

### Control flow inferred from the code

- Line 59: `for` loop over `l=0:5`.
- Line 60: `for` loop over `m=l:-1:-l`.
- Line 78: `for` loop over `l=0:3`.
- Line 79: `for` loop over `m=l:-1:-l`.
- Line 80: `for` loop over `n=l:-1:-l`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/indexing_inverse_suite', 'Indexing inverse helpers', 'indexing helpers must be exact inverses on their finite integer domains.')`.
- Lines 25: computes `S_ref` using `S_ref=[1 3 6 10;2 5 9 13;4 8 12 15;7 11 14 16]`.
- Lines 26: computes `S` using `S=serpentine(4,1)`.
- Lines 36: computes `[K,Q]` using `[K,Q]=ndgrid(1:4,1:4)`.
- Lines 37: computes `I` using `I=kq2lin(4,K,Q,1)`.
- Lines 38: computes `[K_obs,Q_obs]` using `[K_obs,Q_obs]=lin2kq(4,I,1)`.
- Lines 58: computes `L` using `L=[]; M=[]`.
- Lines 61: computes `L(end+1)` using `L(end+1)=l`.
- Lines 62: computes `M(end+1)` using `M(end+1)=m`.
- Lines 68: computes `[L_obs,M_obs]` using `[L_obs,M_obs]=lin2lm(I)`.
- Lines 83: computes `N(end+1)` using `N(end+1)=n`.
- Lines 90: computes `[L_obs,M_obs,N_obs]` using `[L_obs,M_obs,N_obs]=lin2lmn(I)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks serpentine matrix indexing, spin-state L,M indexing,
- and Wigner-function L,M,N indexing over complete low-rank domains.

## Implementation structure

- Tests indexing helper inverses. Syntax:
- result=test_indexing_inverse_suite()
- result -regression test result with explanatory messages
- The test checks serpentine matrix indexing, spin-state L,M indexing,
- and Wigner-function L,M,N indexing over complete low-rank domains.
- Announce the test target
- State the indexing target of the test
- Check the documented base-one serpentine matrix
- Check the documented base-zero serpentine matrix
- Check k,q to linear and back in base-one indexing
- Check k,q to linear and back in base-zero indexing
- Build a complete low-rank L,M domain in documented order

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `serpentine()`, `test_close()`, `kq2lin()`, `lin2kq()`, `lm2lin()`, `lin2lm()`, `lmn2lin()`, `lin2lmn()`.
