# tests/kernel/test_ctx_powder_average.m

- Signature: `result=test_ctx_powder_average()`

## Purpose

Tests powder averaging against explicit weighted summation. Syntax: result=test_ctx_powder_average()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test asks powder() for individual orientation traces and checks that
- the default powder average is the same weighted sum.

## Implementation structure

- Tests powder averaging against explicit weighted summation. Syntax:
- result=test_ctx_powder_average()
- result -regression test result with explanatory messages
- The test asks powder() for individual orientation traces and checks that
- the default powder average is the same weighted sum.
- Announce the test target
- State the powder-averaging target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a tiny powder acquisition
- Run the averaged powder calculation
- Run the per-orientation powder calculation
- Assemble the independent weighted sum
