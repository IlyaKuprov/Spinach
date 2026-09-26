# examples/fundamentals/quadratures/aht_benchmark.m

- Signature: `aht_benchmark()`

## Purpose

Accuracy benchmark for the Hamiltonian period propagator caclulation using Lie group integrators. Bacause the mo- dulation is purely sinusoidal, 2nd and 4th order integra- tors show the same apparent accuracy. Calculation time: seconds

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Accuracy benchmark for the Hamiltonian period propagator
- caclulation using Lie group integrators. Bacause the mo-
- dulation is purely sinusoidal, 2nd and 4th order integra-
- tors show the same apparent accuracy.
- Calculation time: seconds
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Spectrum setup
