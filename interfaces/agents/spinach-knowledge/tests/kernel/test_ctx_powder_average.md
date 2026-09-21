# tests/kernel/test_ctx_powder_average.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_powder_average.m`
- Signature: `result=test_ctx_powder_average()`
- Total lines: 65

## Purpose

Tests powder averaging against explicit weighted summation. Syntax: result=test_ctx_powder_average()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `powder()`, `test_spin_system()`, `state()`, `test_close()`.
