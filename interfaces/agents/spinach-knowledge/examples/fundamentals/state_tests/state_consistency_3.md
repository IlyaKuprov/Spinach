# examples/fundamentals/state_tests/state_consistency_3.m

- Signature: `state_consistency_3()`

## Purpose

Deuterium pair singlet, triplet, and quintet state internal consistency test.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Deuterium pair singlet, triplet, and quintet state
- internal consistency test.
- A pair of deuteria
- Hilbert space
- Spinach housekeeping
- Ortho-deuterium states from Spinach
- Component vectors in Hilbert space as per Eq 1
- in https://doi.org/10.1016/S0009-2614(98)00784-2
- Test singlet state
- Test triplet states
- Test quintet states
- Move to Liouville space
