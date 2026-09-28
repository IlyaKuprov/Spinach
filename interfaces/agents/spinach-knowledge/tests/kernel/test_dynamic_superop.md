# tests/kernel/test_dynamic_superop.m

- Signature: `result=test_dynamic_superop()`

## Purpose

Regression tests for left, right, and commutator superoperators.

## Tests

- Constructs two-spin `[2 2]` product-operator cases and checks `A_comm = A_left - A_right`.
- Checks that the commutator superoperator is nonzero and sparse.
- Checks `local_xyz_to_sparse`, which converts XYZ triples to a MATLAB sparse matrix.

## Outputs

- result -regression test result with explanatory messages
- The test checks the unit-operator shortcut, direct commutator and
- anticommutator identities, and Lz spherical-tensor projection eigenvalues.
