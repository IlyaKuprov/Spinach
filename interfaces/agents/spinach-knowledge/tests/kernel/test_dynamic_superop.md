# tests/kernel/test_dynamic_superop.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_superop.m`
- Signature: `result=test_dynamic_superop()`
- Total lines: 74

## Purpose

Tests superop() sparse-XYZ spherical-tensor product operators. Syntax: result=test_dynamic_superop()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_xyz_to_sparse()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks the unit-operator shortcut, direct commutator and
- anticommutator identities, and Lz spherical-tensor projection eigenvalues.

## Implementation structure

- Tests superop() sparse-XYZ spherical-tensor product operators. Syntax:
- result=test_dynamic_superop()
- result -regression test result with explanatory messages
- The test checks the unit-operator shortcut, direct commutator and
- anticommutator identities, and Lz spherical-tensor projection eigenvalues.
- Announce the test target
- State the superop target of the test
- Build a two-spin spherical-tensor Liouville-space system
- Check the inactive-spin shortcut returns the full identity
- Build left, right, commutator, and anticommutator forms for Lz on spin 1
- Check algebraic side-product identities
- Check the Lz commutator eigenvalues on irreducible tensor projections

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `superop()`, `test_spin_system()`, `local_xyz_to_sparse()`, `test_close()`, `speye()`, `lin2lm()`, `test_true()`, `nnz()`, `xyz()`, `complex()`.
