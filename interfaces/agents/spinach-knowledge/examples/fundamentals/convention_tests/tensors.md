# examples/fundamentals/convention_tests/tensors.m

- Signature: `tensors()`

## Purpose

Test the conversion from Stevens operator coefficients to irreducible spherical tensor coefficients.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

## Implementation structure

- Test the conversion from Stevens operator coefficients to
- irreducible spherical tensor coefficients.
- Test up to rank 6 on spin 15/2
- Generate a random set of Stevens operator coefficients
- Build the linear combination
- Translate the coefficients into ISTs
- Subtract the matrices and check the norm
