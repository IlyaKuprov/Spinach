# tests/kernel/test_indexing_roundtrip_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_indexing_roundtrip_suite.m`
- Signature: `result=test_indexing_roundtrip_suite()`
- Total lines: 64

## Purpose

Tests angular-momentum and matrix indexing helpers. Syntax: result=test_indexing_roundtrip_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

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
