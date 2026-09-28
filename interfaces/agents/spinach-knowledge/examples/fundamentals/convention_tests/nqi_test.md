# examples/fundamentals/convention_tests/nqi_test.m

- Signature: `nqi_test()`

## Purpose

Tests the reverse decomposition of a spin-1 Hamiltonian by reconstructing a random traceless Hermitian 3×3 Hamiltonian from the parameters returned by `ham2nqi`.

## Method and check

The test draws a complex random matrix, Hermitian-symmetrises it, removes its trace, and calls `ham2nqi` to obtain `omega` and `Q`. It then builds a 14N spin system with `Q/(2*pi)` as the coupling matrix, zero magnet field, and the `zeeman-hilb` formalism with no approximation. Spinach reconstructs the Hamiltonian from the quadrupolar term and the three `omega`-weighted `Lx`, `Ly`, and `Lz` operators. The test passes when the 2-norm residual is at most 10⁻⁶ times the 2-norm of the sum of the two Hamiltonians.
