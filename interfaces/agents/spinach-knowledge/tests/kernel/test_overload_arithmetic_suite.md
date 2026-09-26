# tests/kernel/test_overload_arithmetic_suite.m

- Signature: `result=test_overload_arithmetic_suite()`

## Purpose

Tests cheap overload arithmetic for cell, struct, RCV, and polyadic classes. Syntax: result=test_overload_arithmetic_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks elementwise cell overloads, recursive struct arithmetic,
- RCV sparse storage operations, and small polyadic arithmetic against
- explicit Matlab matrix references.

## Implementation structure

- Tests cheap overload arithmetic for cell, struct, RCV, and polyadic classes. Syntax:
- result=test_overload_arithmetic_suite()
- result -regression test result with explanatory messages
- The test checks elementwise cell overloads, recursive struct arithmetic,
- RCV sparse storage operations, and small polyadic arithmetic against
- explicit Matlab matrix references.
- Announce the test target
- State the utility target of the test
- Define small cell-array operands
- Check cell addition and subtraction overloads
- Check cell scalar-array and matrix multiplication overloads
- Check cell totals and inflation shorthand
