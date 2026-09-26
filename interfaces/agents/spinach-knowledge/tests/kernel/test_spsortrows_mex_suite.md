# tests/kernel/test_spsortrows_mex_suite.m

- Signature: `result=test_spsortrows_mex_suite()`

## Purpose

Tests the sparse sortrows MEX helper. Syntax: result=test_spsortrows_mex_suite()

## Physical / mathematical content

## Numerical / algorithmic content

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
