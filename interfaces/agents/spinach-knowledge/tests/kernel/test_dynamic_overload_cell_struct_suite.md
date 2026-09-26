# tests/kernel/test_dynamic_overload_cell_struct_suite.m

- Signature: `result=test_dynamic_overload_cell_struct_suite()`

## Purpose

Tests dynamic dispatch of cheap cell, struct, and double overloads. Syntax: result=test_dynamic_overload_cell_struct_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test exercises object-operation dispatch for cell arithmetic,
- cell multiplication, cell utility overloads, recursive structure
- arithmetic, and the double inflate no-op.

## Implementation structure

- Tests dynamic dispatch of cheap cell, struct, and double overloads. Syntax:
- result=test_dynamic_overload_cell_struct_suite()
- result -regression test result with explanatory messages
- The test exercises object-operation dispatch for cell arithmetic,
- cell multiplication, cell utility overloads, recursive structure
- arithmetic, and the double inflate no-op.
- Announce the test target
- State the overload target of the test
- Build deterministic matrix operands
- Exercise cell plus and minus through operator dispatch
- Exercise numeric-cell plus and minus dispatch
- Exercise direct method-name dispatch for coverage tracking
