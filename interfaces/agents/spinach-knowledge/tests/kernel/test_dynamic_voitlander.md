# tests/kernel/test_dynamic_voitlander.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_voitlander.m`
- Signature: `result=test_dynamic_voitlander()`
- Total lines: 191

## Purpose

Tests voitlander() on an isotropic one-electron field-swept line. Syntax: result=test_dynamic_voitlander()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `voitlander()`, `test_spin_system()`, `state()`, `hamiltonian()`, `assume()`, `orientation()`, `eigenfields()`, `test_true()`, `isscalar()`, `test_close()`, `sphtrsubd()`, `sphtarea()`, `all()`, `spec()`, `tri_corr()`.
