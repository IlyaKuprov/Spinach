# examples/fundamentals/tensor_structures/polyadic_test_1.m

- Signature: `polyadic_test_1()`

## Purpose

Unit tests for the polyadic object.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Unit tests for the polyadic object.
- Get random complex matrices
- Get a random complex vector
- Normalise everything
- Create the polyadic
- Add prefixes and suffixes
- Reference matrix
- Create-inflate test
- Matrix-vector test
- Addition test 1
- Addition test 2
- Conjugate-transpose test 1
