# examples/dnp_liq/jdnp/fig_4_state_amplitudes.m

- Signature: `fig_4_state_amplitudes()`

## Purpose

Time evolution of the individual states in the basis set built from the singlet-triplet basis on the two electrons and Carte- sian spin operator basis on the nucleus, demonstrating the im- balance between singlet-alpha and singlet-beta subspace. This leads to trainsient nuclear spin polarisation. Details in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Time evolution of the individual states in the basis set built
- from the singlet-triplet basis on the two electrons and Carte-
- sian spin operator basis on the nucleus, demonstrating the im-
- balance between singlet-alpha and singlet-beta subspace. This
- leads to trainsient nuclear spin polarisation. Details in:
- Calculation time: seconds
- Load the spin system
- Experiment parameters
- Set the magnet
- Set microwave offset frequency
- Set the exchange coupling
- Spinach housekeeping
