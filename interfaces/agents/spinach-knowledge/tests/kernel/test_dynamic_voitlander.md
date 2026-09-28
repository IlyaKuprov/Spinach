# tests/kernel/test_dynamic_voitlander.m

- Signature: `result=test_dynamic_voitlander()`

## Purpose

Regression tests for the spherical-triangle integration routine.

## Tests

- Defines a spherical triangle and compares the numerical result with an analytic Lorentzian reference.
- Checks vertex-product averaging.
- Checks rejection of an infinite Jacobian and continuity when branches are relabelled.

## Outputs

- result -regression test result with explanatory messages
- The test uses an isotropic spin-half electron, for which all triangle
- vertices and subdivision midpoints have the same transition field.
