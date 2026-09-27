# examples/quantum_tech/tavis_cummings_splitting.m

- Signature: `tavis_cummings_splitting()`

## Purpose

Collective normal-mode splitting in the Tavis-Cummings model for one to four identical electron spins coupled to a common microwave cavity mode. The bright-state splitting follows the square-root scaling of Tavis and Cummings, Phys. Rev. 170, 379 (1968). Calculation time: seconds.

## Model and parameters

The example sets the magnet field to zero and couples one to four identical electron spins (isotope E) to a resonant C3 cavity mode. Each spin-cavity coupling is 3 MHz; the Zeeman-Hilbert basis is used without approximation. The cited Tavis-Cummings model predicts a bright-mode splitting proportional to the square root of the spin count.

## Calculation

For each ensemble size, the code constructs the Hamiltonian, projects it onto the one-excitation manifold, and diagonalizes that subspace. It compares the numerical splitting with `2*sqrt(n)*2*pi*coupling`, failing if the relative discrepancy exceeds 1e-10, then plots both values in MHz. The citation is Tavis and Cummings, Phys. Rev. 170, 379 (1968).
