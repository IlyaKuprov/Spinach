# examples/fundamentals/state_tests/normalization_1.m

- Signature: `normalization_1()`

## Purpose

Internal consistency test for the state vectors and matrices in each of the three formalisms supported by Spinach.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Internal consistency test for the state vectors and matrices
- in each of the three formalisms supported by Spinach.
- System specification
- Preallocate the answer
- Compute norm differences
- Display the answers
