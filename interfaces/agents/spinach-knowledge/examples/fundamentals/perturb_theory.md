# examples/fundamentals/perturb_theory.m

- Signature: `perturb_theory()`

## Purpose

Rayleigh-Schrodinger and Van Vleck perturbation theory modules test. Eigenvector representations differ in the two theories, but the energies are the same.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Rayleigh-Schrodinger and Van Vleck perturbation theory
- modules test. Eigenvector representations differ in the
- two theories, but the energies are the same.
- Settings
- H0 -Zeeman interaction
- H1 -random matrix
- Energies -perturbation theories
- Energies -diagonalisation
- Comparison with diagonalisation
- Eigensystems, PTs
- RSPT gets the eigensystem directly
- VVPT returns a generator that needs exponentiation
