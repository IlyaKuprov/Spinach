# tests/kernel/test_lowlevel_utilities_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_lowlevel_utilities_suite.m`
- Signature: `result=test_lowlevel_utilities_suite()`
- Total lines: 105

## Purpose

Tests cheap deterministic low-level utility functions. Syntax: result=test_lowlevel_utilities_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `comm()`, `rocomm()`, `remtrace()`, `remncomm()`, `hdot()`, `atranspose()`, `killcross()`, `killdiag()`, `keep_rank()`, `frob_chop()`, `gaussfun()`, `lorentzfun()`, `spden()`, `test_true()`.
