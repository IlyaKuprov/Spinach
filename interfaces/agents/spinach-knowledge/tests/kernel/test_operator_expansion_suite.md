# tests/kernel/test_operator_expansion_suite.m

- Signature: `result=test_operator_expansion_suite()`

## Purpose

Tests operator expansion and conversion helpers. Syntax: result=test_operator_expansion_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks that irreducible spherical tensor and bosonic monomial
- expansion helpers reconstruct explicit matrices, and that operator-sized
- allocation helpers return the correct formalism dimensions.

## Implementation structure

- Tests operator expansion and conversion helpers. Syntax:
- result=test_operator_expansion_suite()
- result -regression test result with explanatory messages
- The test checks that irreducible spherical tensor and bosonic monomial
- expansion helpers reconstruct explicit matrices, and that operator-sized
- allocation helpers return the correct formalism dimensions.
- Announce the test target
- State the expansion target of the test
- Check Hilbert-to-Liouville vectorisation identities on a non-diagonal matrix
- Check IST expansion of a generic spin-one matrix
- Check spin and boson energy-level counting conventions in IST expansions
- Check central-transition and boson-product IST expansion wrappers
