# tests/kernel/test_indexing_roundtrip_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_indexing_roundtrip_suite.m`
- Signature: `result=test_indexing_roundtrip_suite()`
- Total lines: 64

## Purpose

Tests angular-momentum and matrix indexing helpers. Syntax: result=test_indexing_roundtrip_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Indexing conversion functions\n')`.
- Lines 20-23: State the indexing target of the test; implemented by `result=new_test_result('kernel/indexing_roundtrip_suite', 'Indexing conversion functions', 'indexing helpers must be exact inverses on valid integer domains.')`.
- Lines 25-26: Linear spin-state indexing is zero-based and ordered by increasing L; implemented by `I=0:24`.
- Lines 35-36: Wigner D-function indexing is one-based and ordered by increasing L, then M, then N; implemented by `J=1:35`.
- Lines 44-45: Serpentine matrix indexing has documented triangular scan order; implemented by `S1=serpentine(4,1)`.
- Lines 53-54: Serpentine k,q coordinates and linear indices must be exact inverses in both bases; implemented by `N=4`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/indexing_roundtrip_suite', 'Indexing conversion functions', 'indexing helpers must be exact inverses on valid integer domains.')`.
- Lines 26: computes `I` using `I=0:24`.
- Lines 27: computes `[L,M]` using `[L,M]=lin2lm(I)`.
- Lines 36: computes `J` using `J=1:35`.
- Lines 37: computes `[Lw,Mw,Nw]` using `[Lw,Mw,Nw]=lin2lmn(J)`.
- Lines 45: computes `S1` using `S1=serpentine(4,1)`.
- Lines 46: computes `S0` using `S0=serpentine(4,0)`.
- Lines 47: computes `S1_ref` using `S1_ref=[1 3 6 10; 2 5 9 13; 4 8 12 15; 7 11 14 16]`.
- Lines 54: computes `N` using `N=4`.
- Lines 55: computes `idx1` using `idx1=1:N^2`.
- Lines 56: computes `[K1,Q1]` using `[K1,Q1]=lin2kq(N,idx1,1)`.
- Lines 59: computes `idx0` using `idx0=0:(N^2-1)`.
- Lines 60: computes `[K0,Q0]` using `[K0,Q0]=lin2kq(N,idx0,0)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks that linear and structured index representations are
- mutually consistent for spherical tensors, Wigner functions, and matrix
- serpentine indexing.

## Implementation structure

- Tests angular-momentum and matrix indexing helpers. Syntax:
- result=test_indexing_roundtrip_suite()
- result -regression test result with explanatory messages
- The test checks that linear and structured index representations are
- mutually consistent for spherical tensors, Wigner functions, and matrix
- serpentine indexing.
- Announce the test target
- State the indexing target of the test
- Linear spin-state indexing is zero-based and ordered by increasing L
- Wigner D-function indexing is one-based and ordered by increasing L, then M, then N
- Serpentine matrix indexing has documented triangular scan order
- Serpentine k,q coordinates and linear indices must be exact inverses in both bases

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `lin2lm()`, `test_close()`, `lm2lin()`, `lin2lmn()`, `lmn2lin()`, `serpentine()`, `lin2kq()`, `kq2lin()`.
