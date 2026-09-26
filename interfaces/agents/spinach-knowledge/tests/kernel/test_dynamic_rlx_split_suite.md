# tests/kernel/test_dynamic_rlx_split_suite.m

- Signature: `result=test_dynamic_rlx_split_suite()`

## Purpose

Tests relaxation-superoperator component splitting. Syntax: result=test_dynamic_rlx_split_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks that rlx_split() separates single-spin longitudinal,
- single-spin transverse, and multi-spin relaxation blocks without overlap.

## Implementation structure

- Tests relaxation-superoperator component splitting. Syntax:
- result=test_dynamic_rlx_split_suite()
- result -regression test result with explanatory messages
- The test checks that rlx_split() separates single-spin longitudinal,
- single-spin transverse, and multi-spin relaxation blocks without overlap.
- Announce the test target
- State the relaxation-splitting target of the test
- Build a two-spin spherical-tensor Liouville basis
- Interpret the basis using the same documented state categories
- Build a diagonal relaxation matrix with a zero unit-state element
- Build the mathematically expected non-overlapping blocks
- Split the relaxation superoperator with the production helper
