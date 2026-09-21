# tests/kernel/test_states_composite_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_states_composite_suite.m`
- Signature: `result=test_states_composite_suite()`
- Total lines: 125

## Purpose

Tests composite state generators in kernel/states. Syntax: result=test_states_composite_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks unit-state normalisation, two-spin singlet/triplet
- projectors, four-spin product states, and partner-state enumeration.

## Implementation structure

- Tests composite state generators in kernel/states. Syntax:
- result=test_states_composite_suite()
- result -regression test result with explanatory messages
- The test checks unit-state normalisation, two-spin singlet/triplet
- projectors, four-spin product states, and partner-state enumeration.
- Announce the test target
- State the state-generation target of the test
- Build a one-proton Hilbert-space spin system
- Check the Hilbert-space thermodynamic unit state
- Check the Zeeman-Liouville unit vector normalisation
- Build a two-proton Hilbert-space spin system
- Textbook two-spin wavefunctions in the Zeeman product basis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `test_close()`, `unit_state()`, `speye()`, `unit()`, `singlet()`, `triplet()`, `four_spin_states()`, `partner_state()`, `test_true()`, `isequal()`, `descr()`, `expected_descr()`, `state()`, `int2str()`.
