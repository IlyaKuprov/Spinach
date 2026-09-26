# tests/kernel/test_commutator_utility.m

- Signature: `result=test_commutator_utility()`

## Purpose

Tests the commutator utility. Syntax: result=test_commutator_utility()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks that comm(A,B) implements AB-BA and returns zero for
- mutually commuting matrices.

## Implementation structure

- Tests the commutator utility. Syntax:
- result=test_commutator_utility()
- result -regression test result with explanatory messages
- The test checks that comm(A,B) implements AB-BA and returns zero for
- mutually commuting matrices.
- Announce the test target
- State the mathematical target of the test
- Define a non-commuting pair and its reference commutator
- Check the utility on a non-commuting pair
- Check a commuting diagonal pair
