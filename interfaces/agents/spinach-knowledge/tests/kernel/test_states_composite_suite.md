# tests/kernel/test_states_composite_suite.m

- Signature: `result=test_states_composite_suite()`

## Purpose

Tests composite state generators in kernel/states. Syntax: result=test_states_composite_suite()

## Physical / mathematical content

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
