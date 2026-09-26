# tests/kernel/test_indexing_inverse_suite.m

- Signature: `result=test_indexing_inverse_suite()`

## Purpose

Tests indexing helper inverses. Syntax: result=test_indexing_inverse_suite()

## Physical / mathematical content

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
