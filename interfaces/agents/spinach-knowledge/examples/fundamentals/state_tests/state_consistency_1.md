# examples/fundamentals/state_tests/state_consistency_1.m

- Signature: `state_consistency_1()`

## Purpose

Test of internal consistency for state and operator generation across the three formalisms supported by Spinach. Two-and four-spin states are tested.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Test of internal consistency for state and operator
- generation across the three formalisms supported by
- Spinach. Two-and four-spin states are tested.
- Magnet field
- Set the spin system
- Loop over the formalisms
- Basis set
- Hush the logs
- Spinach housekeeping
- Unit state from Spinach
- Two-spin singlet-triplet state sum test
- Four-spin singlet-triplet state sum test
