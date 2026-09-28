# tests/kernel/test_perm_group_suite.m

- Signature: `result=test_perm_group_suite()`

## Purpose

Tests permutation group database metadata.

## Physical / mathematical content

- Tests real-valued-character Abelian permutation subgroups, including character orthogonality, closure, commutativity, and involution order.

## Numerical / algorithmic content

- Compares the S6A and S8A element and character tables with explicit reference tables using tolerances of `1e-15`.
- Checks S4A, S6A, and S8A metadata, singleton classes, one-dimensional irreducible representations, real sign-valued characters, valid permutation rows, and maximality via centraliser size.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announces the test target and initialises the regression result.
- Checks the S6A and S8A element and character tables.
- Checks structural consistency of S4A, S6A, and S8A, including permutation products and centralisers.