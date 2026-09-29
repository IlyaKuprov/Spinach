# tests/kernel/test_commutator_utility.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_commutator_utility.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_commutator_utility.m)

## Purpose

Regression test for the matrix commutator utility `comm(A,B)`. The test verifies that `comm(A,B)` implements the definition AB−BA exactly, and that it returns zero for mutually commuting matrices.

## Behaviour

- Announces the test target by printing `TESTING: Matrix commutator utility`.
- Initialises a regression test result via `new_test_result` with suite `kernel/commutator_utility`, name `Matrix commutator utility`, and the target statement that `comm(A,B)` must return AB−BA exactly.
- Defines a non-commuting pair `A=[1 2;3 4]` and `B=[0 1;-1 2]`, computes the reference commutator `C=A*B-B*A`, and checks `comm(A,B)` against `C` using `test_close` with absolute and relative tolerances of `1e-15`, under the label `non-commuting pair` with the message `the commutator is defined as AB-BA`.
- Defines a commuting diagonal pair `D=diag([1 2 3])` and `E=diag([4 5 6])`, and checks `comm(D,E)` against `zeros(3)` using `test_close` with absolute and relative tolerances of `1e-15`, under the label `commuting diagonal pair` with the message `diagonal matrices in the same basis commute`.

## Inputs and outputs

- Inputs: none. The function is called as `result=test_commutator_utility()`.
- Outputs:
  - `result` — regression test result with explanatory messages, accumulated through successive `test_close` checks.

## References

- `comm` — matrix commutator utility under test.
- `new_test_result` — initialises the regression test result.
- `test_close` — performs the closeness checks against reference values.

Contact: ilya.kuprov@weizmann.ac.il (per the source header).
