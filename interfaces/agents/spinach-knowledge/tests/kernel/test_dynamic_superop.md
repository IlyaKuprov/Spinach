# tests/kernel/test_dynamic_superop.m

- Signature: `result=test_dynamic_superop()`

## Purpose

Tests superop() sparse-XYZ spherical-tensor product operators. Syntax: result=test_dynamic_superop()

## Physical / mathematical content

## Numerical / algorithmic content

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
