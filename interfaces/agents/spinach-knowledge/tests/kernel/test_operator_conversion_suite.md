# tests/kernel/test_operator_conversion_suite.m

- Signature: `result=test_operator_conversion_suite()`

## Purpose

Tests Hilbert-to-Liouville operator conversion utilities. Syntax: result=test_operator_conversion_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks direct vectorisation identities for left, right,
- commutation, and anticommutation superoperators, unit_oper dimensions,
- and Lindbladian rate calibration.

## Implementation structure

- Tests Hilbert-to-Liouville operator conversion utilities. Syntax:
- result=test_operator_conversion_suite()
- result -regression test result with explanatory messages
- The test checks direct vectorisation identities for left, right,
- commutation, and anticommutation superoperators, unit_oper dimensions,
- and Lindbladian rate calibration.
- Announce the test target
- State the operator-conversion target of the test
- Define a non-trivial Hermitian operator
- Check direct vectorisation formulas
- Check unit operator dimensions in major formalisms
- Check Lindbladian calibration on a simple vector
