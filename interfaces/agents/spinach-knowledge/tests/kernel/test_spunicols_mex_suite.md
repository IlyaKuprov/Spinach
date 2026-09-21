# tests/kernel/test_spunicols_mex_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_spunicols_mex_suite.m`
- Signature: `result=test_spunicols_mex_suite()`
- Total lines: 92

## Purpose

Tests the sparse unique-column MEX helper. Syntax: result=test_spunicols_mex_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

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
