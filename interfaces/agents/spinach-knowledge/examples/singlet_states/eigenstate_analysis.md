# examples/singlet_states/eigenstate_analysis.m

- Signature: `eigenstate_analysis()`

## Purpose

Stationary state analysis for the spin system of allyl pyruvate, finding out which component of the singlet state commutes with the drift Hamiltonian.

## Physical / mathematical content

For allyl pyruvate at 14.1 T, the analysis tests which part of the singlet state on spins 3 and 4 commutes with the drift Hamiltonian, then repeats with a 2850 offset and a 1 kHz 1H spin-lock.

## Numerical / algorithmic content

Each pass diagonalizes the Hermitian Hamiltonian, removes noncommuting and unit components from the singlet state, projects out the normalized Lz–Lz component, and reports the remaining norms.

## Implementation structure

The code builds a 1H/13C spin system from `allyl_pyruvate`, selects the fourth single-13C isotopomer, uses an unapproximated `zeeman-hilb` basis and isotropic NMR Hamiltonian, then repeats the analysis after a 2850 offset and a `2*pi*1000` 1H `Lx` term.
