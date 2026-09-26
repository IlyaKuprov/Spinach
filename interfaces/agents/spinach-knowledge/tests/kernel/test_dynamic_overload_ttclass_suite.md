# tests/kernel/test_dynamic_overload_ttclass_suite.m

- Signature: `result=test_dynamic_overload_ttclass_suite()`

## Purpose

Tests dynamic dispatch of tensor-train class overloads. Syntax: result=test_dynamic_overload_ttclass_suite()

## Physical / mathematical content

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Outputs

- result -regression test result with explanatory messages
- The test exercises ttclass object operations against exact dense
- references on one-core and two-core tensor trains.

## Implementation structure

- Tests dynamic dispatch of tensor-train class overloads. Syntax:
- result=test_dynamic_overload_ttclass_suite()
- result -regression test result with explanatory messages
- The test exercises ttclass object operations against exact dense
- references on one-core and two-core tensor trains.
- Announce the test target
- State the overload target of the test
- Build one-core tensor trains and dense references
- Exercise constructor, full, size, sizes, ranks, numel, and subsref dispatch
- Exercise addition, subtraction, scalar multiplication, and scalar division
- Exercise matrix, vector, and tensor-train multiplication dispatch
- Exercise conjugation, transposition, trace, diagonal, sum, and mean dispatch
