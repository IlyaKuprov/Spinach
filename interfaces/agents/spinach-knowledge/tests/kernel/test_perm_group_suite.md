# tests/kernel/test_perm_group_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_perm_group_suite.m`
- Signature: `result=test_perm_group_suite()`
- Total lines: 142

## Purpose

Tests permutation group database metadata. Syntax: result=test_perm_group_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `perm_group()`, `test_close()`, `test_true()`, `group_orders()`, `group_degrees()`, `all()`, `isequal()`, `element()`, `ismember()`, `perms()`, `all_perms()`.
