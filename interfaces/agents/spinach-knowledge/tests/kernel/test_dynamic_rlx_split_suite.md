# tests/kernel/test_dynamic_rlx_split_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_rlx_split_suite.m`
- Signature: `result=test_dynamic_rlx_split_suite()`
- Total lines: 69

## Purpose

Tests relaxation-superoperator component splitting. Syntax: result=test_dynamic_rlx_split_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `rlx_split()`, `test_spin_system()`, `lin2lm()`, `logical()`, `any()`, `diag_vals()`, `spdiags()`, `R1_ref()`, `R2_ref()`, `Rm_ref()`, `test_close()`.
