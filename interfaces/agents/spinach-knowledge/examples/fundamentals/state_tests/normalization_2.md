# examples/fundamentals/state_tests/normalization_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/state_tests/normalization_2.m`
- Signature: `normalization_2()`
- Total lines: 52

## Purpose

Internal consistency test for the state vectors and matrices. Checks that the inner products are consistent between the for- malisms supported by Spinach. Output should be a 6x3 matrix with identical columns.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Internal consistency test for the state vectors and matrices.
- Checks that the inner products are consistent between the for-
- malisms supported by Spinach. Output should be a 6x3 matrix
- with identical columns.
- System specification
- Preallocate the answer
- Get the norms
- Run the tests

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `norms()`.
