# tests/kernel/test_indexing_inverse_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_indexing_inverse_suite.m`
- Signature: `result=test_indexing_inverse_suite()`
- Total lines: 102

## Purpose

Tests indexing helper inverses. Syntax: result=test_indexing_inverse_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

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
