# tests/kernel/test_lowlevel_utilities_suite.m

- Signature: `result=test_lowlevel_utilities_suite()`

## Purpose

Tests cheap deterministic low-level utility functions. Syntax: result=test_lowlevel_utilities_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks small numerical helpers, matrix filters, integer type
- selection, and analytic line-shape definitions against explicit answers.

## Implementation structure

- Tests cheap deterministic low-level utility functions. Syntax:
- result=test_lowlevel_utilities_suite()
- result -regression test result with explanatory messages
- The test checks small numerical helpers, matrix filters, integer type
- selection, and analytic line-shape definitions against explicit answers.
- Announce the test target
- State the utility target of the test
- Define small test matrices
- Check commutator and right-ordered nested commutator
- Check trace removal and commuting part extraction
- Check Frobenius inner product and anti-diagonal transpose
- Check matrix wiping helpers
