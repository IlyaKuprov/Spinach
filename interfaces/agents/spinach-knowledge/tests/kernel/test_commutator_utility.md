# tests/kernel/test_commutator_utility.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_commutator_utility.m`
- Signature: `result=test_commutator_utility()`
- Total lines: 40

## Purpose

Tests the commutator utility. Syntax: result=test_commutator_utility()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `comm()`, `test_close()`.
