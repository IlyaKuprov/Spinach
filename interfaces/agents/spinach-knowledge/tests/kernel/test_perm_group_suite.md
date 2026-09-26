# tests/kernel/test_perm_group_suite.m

- Signature: `result=test_perm_group_suite()`

## Purpose

Tests permutation group database metadata. Syntax: result=test_perm_group_suite()

## Physical / mathematical content

- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks the real-valued-character Abelian permutation subgroups
- against explicit element tables, explicit character tables, character
- orthogonality, closure, commutativity, and involution order.

## Implementation structure

- Tests permutation group database metadata. Syntax:
- result=test_perm_group_suite()
- result -regression test result with explanatory messages
- The test checks the real-valued-character Abelian permutation subgroups
- against explicit element tables, explicit character tables, character
- orthogonality, closure, commutativity, and involution order.
- Announce the test target
- State the utility target of the test
- Check the S6 Abelian subgroup table
- Check the S8 Abelian subgroup table
- Check structural consistency of all real-valued-character Abelian options
- Get the group under test
