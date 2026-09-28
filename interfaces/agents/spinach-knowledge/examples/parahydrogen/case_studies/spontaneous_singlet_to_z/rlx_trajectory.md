# examples/parahydrogen/case_studies/spontaneous_singlet_to_z/rlx_trajectory.m

- Signature: `rlx_trajectory()`

## Purpose

Tracks the time dependence of three spin orders in a para-hydrogen molecule coordinated to a nickel cage with large chemical-shift anisotropy. The source notes that the paper link is forthcoming and points to a Mathematica worksheet.

## Physical / mathematical content

The two-spin system contains `1H` nuclei with supplied coordinates and trace-subtracted anisotropic Zeeman tensors. It starts in the two-spin singlet state and evolves under secular Redfield relaxation, with zero equilibrium state and a 500 ps correlation time. The three detected channels are the transverse exchange order `A/2 = (L+S- + L-S+)/2`, the longitudinal product `4 LzSz`, and total longitudinal order `Lz + Sz`.

## Numerical / algorithmic content

The Liouville-space generator combines the Hamiltonian and relaxation superoperators. A multichannel propagation computes the three observables over 2 ms using 1000 time points.

## Implementation structure

The script sets the field to 18.7893 T, defines both proton coordinates and Zeeman tensors, selects the spherical-tensor Liouville formalism with no basis approximation, constructs the spin system, and then builds the Liouvillian, singlet initial state, detection operators, trajectory, and plot.
