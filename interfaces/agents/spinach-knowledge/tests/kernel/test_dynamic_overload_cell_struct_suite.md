# tests/kernel/test_dynamic_overload_cell_struct_suite.m

## Purpose

Regression test for dynamic dispatch of the cheap cell, struct, and double overloads in Spinach. The suite verifies that operator and function dispatch on small deterministic objects matches explicit dense references, covering cell arithmetic, cell multiplication, cell utility overloads, recursive structure arithmetic, and the double inflate no-op.

Source: [tests/kernel/test_dynamic_overload_cell_struct_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_overload_cell_struct_suite.m)

## Behaviour

- Announces the test target with `fprintf('TESTING: Cell and structure overload dispatch\n')` and initialises a test result via `new_test_result` for `kernel/dynamic_overload_cell_struct_suite`, describing the requirement that cell, struct, and double overloads match explicit dense references on small deterministic objects.
- Builds deterministic 2x2 matrix operands `A=[1 2;3 4]`, `B=[0 5;-1 2]`, `C=[2 -3;4 1]`, `D=[-2 0;5 3]`, and cells `cell_a={A,B}`, `cell_b={C,D}`.
- Cell plus and minus: `cell_a+cell_b` and `cell_a-cell_b` are checked element-wise against `A+C`, `B+D`, `A-C`, `B-D` with tolerances `1e-15`.
- Numeric-cell dispatch: `10+cell_a`, `cell_a-3`, and `7-cell_a` are checked against `10+A`, `B-3`, and `7-A` respectively.
- Direct method-name dispatch: `plus(cell_a,cell_b)` and `minus(cell_a,cell_b)` are checked for coverage tracking of the cell overloads.
- Cell times and mtimes: with `weights=[2 3]` and `R=[2 1;0 -1]`, checks `cell_a.*weights`, `weights.*cell_a`, `cell_a*R`, `R*cell_a`, `times(cell_a,weights)`, and `mtimes(cell_a,R)` against `3*B`, `2*A`, `A*R`, `R*B`, `2*A`, and `B*R`.
- Cell utility overloads: `totsum({sparse(A),sparse(B)})` against `sparse(A+B)`; `complex(cell_a)` against `complex(A)`; `inflate(cell_a)` against `B` (numeric cell entries unchanged); `blkdiag({A},{B})` checked for upper block `A`, lower block `B`, and empty off-diagonal cells, with an error raised if off-diagonal cells are non-empty and a pass message appended otherwise.
- Recursive struct arithmetic: with nested structs `struct_a` and `struct_b` (fields `alpha`, `beta.gamma`, `beta.delta={A,B}` / `{C,D}`), checks `struct_a+struct_b` against `[6;8]`, `[10 12]`, and `A+C`; `2*struct_a` against `2*struct_a.alpha` and `2*B`; and direct `plus(struct_a,struct_b)` and `mtimes(2,struct_a)` dispatch.
- Double inflate no-op: `inflate(A)` is checked against `A` with tolerances `1e-15`, verifying dense numeric arrays are unchanged.
- All comparisons use `test_close` with absolute and relative tolerances of `1e-15` and explanatory messages.

## Inputs and outputs

- **result** - regression test result with explanatory messages, returned by the function.

## References

- [Spinach source: tests/kernel/test_dynamic_overload_cell_struct_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_overload_cell_struct_suite.m)
