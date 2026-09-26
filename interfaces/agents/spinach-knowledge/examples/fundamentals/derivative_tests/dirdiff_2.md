# examples/fundamentals/derivative_tests/dirdiff_2.m

- Signature: `dirdiff_2()`

## Purpose

Test of matrix exponential differentiation of second order Magnus product quadrature (trapdiff.m) with the result com- pared to the central finite difference derivative. General coherent + non-symmetric dissipative case is tested.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Test of matrix exponential differentiation of second order
- Magnus product quadrature (trapdiff.m) with the result com-
- pared to the central finite difference derivative. General
- coherent + non-symmetric dissipative case is tested.
- Formalisms to test
- Loop over formalisms
- Get the Spinach object
- Left and right drift generators, dissipative
- Control operator
- A reasonable time step estimate
- Reasonable controls
- Get analytical derivatives
