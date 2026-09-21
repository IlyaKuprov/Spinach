# tests/kernel/test_state_constructor_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_state_constructor_suite.m`
- Signature: `result=test_state_constructor_suite()`
- Total lines: 120

## Purpose

Tests state-constructor helper functions. Syntax: result=test_state_constructor_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks unit, thermal, singlet/triplet, partner, deuteron-pair,
- four-spin, and zero-field triplet state constructors using projector and
- normalisation identities.

## Implementation structure

- Tests state-constructor helper functions. Syntax:
- result=test_state_constructor_suite()
- result -regression test result with explanatory messages
- The test checks unit, thermal, singlet/triplet, partner, deuteron-pair,
- four-spin, and zero-field triplet state constructors using projector and
- normalisation identities.
- Announce the test target
- State the state-constructor target of the test
- One-spin Hilbert and Liouville unit states have known forms
- Two-spin singlet and triplet constructors must form four orthogonal projectors summing to identity
- partner_state must enumerate all requested partner-state combinations
- Four-spin singlet-singlet state must match the explicit product of two two-spin singlets

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `test_close()`, `unit_state()`, `speye()`, `equilibrium()`, `unit_zeeman()`, `stateinfo()`, `test_true()`, `singlet()`, `triplet()`, `partner_state()`, `int2str()`, `state()`, `four_spin_states()`, `deut_pair()`.
