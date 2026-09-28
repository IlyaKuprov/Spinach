# tests/kernel/test_dynamic_overload_cell_struct_suite.m

- Signature: `result=test_dynamic_overload_cell_struct_suite()`

## Purpose

Regression test for dynamic dispatch of cell, struct, and double overloads against explicit references built from small deterministic objects. The result contains explanatory messages.

## Tests

- Checks cell addition and subtraction, including numeric–cell operations and direct `plus()` and `minus()` calls.
- Checks cell element-wise scaling and left and right matrix multiplication, including direct `times()` and `mtimes()` calls.
- Checks cell utilities: `totsum` of sparse entries, `complex`, `inflate` of numeric entries, and `blkdiag` with preserved diagonal blocks and empty off-diagonal cells.
- Checks recursive struct addition and scalar multiplication across top-level numeric fields, nested numeric fields, and nested cells, including direct method-name dispatch.
- Checks that `inflate` leaves a dense double array unchanged.

Numeric comparisons use absolute and relative tolerances of `1e-15`.