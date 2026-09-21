# tests/kernel/test_spsortrows_mex_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_spsortrows_mex_suite.m`
- Signature: `result=test_spsortrows_mex_suite()`
- Total lines: 84

## Purpose

Tests the sparse sortrows MEX helper. Syntax: result=test_spsortrows_mex_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `sortrows()`, `test_true()`, `isequal()`, `spsortrows()`, `isequaln()`, `onCleanup()`, `rng()`, `randi()`, `sprandn()`, `clear()`.
