# tests/kernel/test_dynamic_cubic_mex_suite.m

- Signature: `result=test_dynamic_cubic_mex_suite()`

## Purpose

Tests the cubic-polynomial MEX helper used by eigenfields(). Syntax: result=test_dynamic_cubic_mex_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks cubic roots, degenerate lower-order polynomials,
- repeated roots, endpoint roots, extreme coefficient scaling, and
- derivative-root use cases against explicit references.

## Implementation structure

- Tests the cubic-polynomial MEX helper used by eigenfields(). Syntax:
- result=test_dynamic_cubic_mex_suite()
- result -regression test result with explanatory messages
- The test checks cubic roots, degenerate lower-order polynomials,
- repeated roots, endpoint roots, extreme coefficient scaling, and
- derivative-root use cases against explicit references.
- Announce the test target
- State the utility target of the test
- Set the production root tolerance
- Check three roots, including endpoints
- Check a triple root
- Check a double root plus an endpoint root
