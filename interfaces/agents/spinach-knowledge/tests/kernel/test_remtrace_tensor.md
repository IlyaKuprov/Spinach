# tests/kernel/test_remtrace_tensor.m

- Signature: `result=test_remtrace_tensor()`

## Purpose

Tests removal of the isotropic tensor trace. Syntax: result=test_remtrace_tensor()

## Physical / mathematical content

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
