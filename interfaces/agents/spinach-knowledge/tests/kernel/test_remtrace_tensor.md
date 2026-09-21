# tests/kernel/test_remtrace_tensor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_remtrace_tensor.m`
- Signature: `result=test_remtrace_tensor()`
- Total lines: 38

## Purpose

Tests removal of the isotropic tensor trace. Syntax: result=test_remtrace_tensor()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks that remtrace subtracts the isotropic component of a
- second-rank interaction tensor, leaving the anisotropic traceless part.

## Implementation structure

- Tests removal of the isotropic tensor trace. Syntax:
- result=test_remtrace_tensor()
- result -regression test result with explanatory messages
- The test checks that remtrace subtracts the isotropic component of a
- second-rank interaction tensor, leaving the anisotropic traceless part.
- Announce the test target
- State the physical target of the test
- Define a symmetric interaction tensor with non-zero isotropic part
- Check explicit trace removal and invariants

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `remtrace()`, `test_close()`, `A_obs()`.
