# tests/kernel/test_indexing_roundtrip_suite.m

- Signature: `result=test_indexing_roundtrip_suite()`

## Purpose

Tests angular-momentum and matrix indexing helpers. Syntax: result=test_indexing_roundtrip_suite()

## Physical / mathematical content

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
