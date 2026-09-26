# examples/fundamentals/derivative_tests/dirdiff_1.m

- Signature: `dirdiff_1()`

## Purpose

Test of matrix exponential differentiation routines. Analytical derivatives are compared to central finite differences.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Test of matrix exponential differentiation routines. Analytical
- derivatives are compared to central finite differences.
- Formalisms to test
- Loop over formalisms
- Get the Spinach object
- Random Hamiltonian
- Random direction operators
- First derivative, numerical
- First derivative, analytical
- Test the first derivative
- Second derivative, numerical
- Second derivative, analytical
