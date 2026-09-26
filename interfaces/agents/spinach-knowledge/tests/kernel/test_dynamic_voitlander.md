# tests/kernel/test_dynamic_voitlander.m

- Signature: `result=test_dynamic_voitlander()`

## Purpose

Tests voitlander() on an isotropic one-electron field-swept line. Syntax: result=test_dynamic_voitlander()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test uses an isotropic spin-half electron, for which all triangle
- vertices and subdivision midpoints have the same transition field.

## Implementation structure

- Tests voitlander() on an isotropic one-electron field-swept line. Syntax:
- result=test_dynamic_voitlander()
- result -regression test result with explanatory messages
- The test uses an isotropic spin-half electron, for which all triangle
- vertices and subdivision midpoints have the same transition field.
- Announce the test target
- State the Voitlander target of the test
- Build an isotropic one-electron Hilbert-space system
- Set field-swept EPR parameters with a deliberately loose recursion gate
- Get Zeeman, coupling, and microwave Hamiltonians
- Find the isotropic transition at one orientation
- Define the positive octant spherical triangle
