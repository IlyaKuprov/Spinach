# examples/fundamentals/perturb_theory.m

- Signature: `perturb_theory()`

## Purpose

Compare Rayleigh-Schrödinger (RSPT) and Van Vleck (VVPT) perturbation energies and eigensystems with direct diagonalisation, through perturbation order 10.

## Physical / mathematical content

A 512-dimensional Zeeman Hamiltonian is perturbed by a random Hermitian matrix scaled by 1/25. The example compares the perturbative energy and eigenvector representations of the two theories with the exact eigensystem of the perturbed Hamiltonian.

## Numerical / algorithmic content

At orders 1:10, energy errors are normalized by the exact energy-vector norm; eigenvector agreement is measured by the 2-norm of abs(V'*V_inf)-I. Both residuals are plotted against perturbation order on logarithmic vertical axes.

## Implementation structure

- Set H0 from the z component of pauli(512) and construct the scaled random Hermitian perturbation H1.
- Call rspert and vvpert at each order, and compare their energies against the sorted exact eigenvalues of H0+H1.
- Obtain the RSPT eigensystem directly; exponentiate the VVPT generator before comparing both eigensystems with exact eigenvectors.
