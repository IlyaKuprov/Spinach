# tests/kernel/test_perm_group_suite.m

## Purpose

Regression test for the permutation group metadata database. It verifies that the real-valued-character Abelian permutation subgroups returned by `perm_group` are exact and internally consistent.

## Behaviour

- Announces the test target with `fprintf` and initialises the result via `new_test_result` under the suite name `kernel/perm_group_suite`, describing the requirement that real-valued-character Abelian permutation subgroup tables be exact and internally consistent.
- Checks the `S6A` table against an explicit 8-element reference table (permutations of degree 6) and an explicit 8-by-8 reference character table, using `test_close` with absolute and relative tolerances of `1e-15`. The reference notes state that `S6A` is represented by the direct product of `S4A` and `S2`, and has the real character table of `C2 x C2 x C2`.
- Checks the `S8A` table against an explicit 16-element reference table (permutations of degree 8) and an explicit 16-by-16 reference character table, using `test_close` with tolerances `1e-15`. The reference notes state that `S8A` is represented by a maximal `C2 x C2 x C2 x C2` subgroup with the corresponding real character table.
- Loops over the group names `{'S4A','S6A','S8A'}` with degrees `[4 6 8]` and orders `[4 8 16]`, and for each group:
  - Verifies scalar metadata: `G.order`, `G.nclasses`, and `G.n_irreps` all equal the group order, and `G.order == 2^floor(degree/2)` (maximal real-character order for the permutation degree).
  - Verifies with `test_close` (tolerances `1e-15`) that `G.class_sizes` and `G.irrep_dims` are all ones (singleton conjugacy classes and one-dimensional irreducible representations).
  - Verifies that `G.class_characters` is real with all entries of absolute value 1 (characters are signs), and that `G.class_characters * G.class_characters.'` equals `order * eye(order)` (row orthogonality).
  - Checks every element row is a valid permutation of `1:degree`, is an involution (`element(element)` equals the identity), and that all pairwise compositions are closed within the listed rows and commute.
  - Enumerates all permutations of `1:degree` with `perms` and counts those commuting with every group element, requiring the centraliser count to equal the group order (the centraliser is the subgroup itself).
- Each check appends a pass/fail entry with an explanatory message to the returned result.

## Inputs and outputs

- **Outputs**:
  - `result` - regression test result object with explanatory messages, produced by `new_test_result` and accumulated through `test_close` and `test_true`.
- The function takes no inputs.

## References

- Source: [tests/kernel/test_perm_group_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_perm_group_suite.m)
- Uses Spinach functions: `new_test_result`, `test_close`, `test_true`, `perm_group`.
