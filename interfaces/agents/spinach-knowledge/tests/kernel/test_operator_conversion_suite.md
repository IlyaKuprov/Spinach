# tests/kernel/test_operator_conversion_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_operator_conversion_suite.m`
- Signature: `result=test_operator_conversion_suite()`
- Total lines: 67

## Purpose

Tests Hilbert-to-Liouville operator conversion utilities. Syntax: result=test_operator_conversion_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `pauli()`, `speye()`, `test_close()`, `hilb2liouv()`, `transpose()`, `test_spin_system()`, `unit_oper()`, `lindbladian()`.
