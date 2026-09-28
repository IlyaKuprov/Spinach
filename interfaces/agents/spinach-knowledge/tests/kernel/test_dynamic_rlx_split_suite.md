# tests/kernel/test_dynamic_rlx_split_suite.m

- Signature: `result=test_dynamic_rlx_split_suite()`

## Purpose

Tests relaxation-superoperator component splitting. Syntax: result=test_dynamic_rlx_split_suite()



## Outputs

- `result` — regression test result with explanatory messages.
## Implementation structure

- Build a two-spin spherical-tensor Liouville basis and classify the single-spin longitudinal, single-spin transverse, and multi-spin components.
- Construct a diagonal relaxation matrix with its unit-state element set to zero, then build the expected longitudinal, transverse, and mixed blocks from the state-category masks.
- Call `rlx_split()`, compare each returned block with its reference, and verify that the blocks reconstruct the original matrix.