# examples/dnp_liq/jdnp/fig_6_microwave_free.m

- Signature: `fig_6_microwave_free()`

## Purpose

A demonstration of Maria Grazia Concilio's microwave-free JDNP effect where a field ramp in combination with unequal relaxati- on rates of singlet-alpha and singlet-beta product states crea- tes nuclear magnetisation enhancement beyond the Boltzmann le- vel at both the starting and the final field. Calculation time: minutes

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A demonstration of Maria Grazia Concilio's microwave-free JDNP
- effect where a field ramp in combination with unequal relaxati-
- on rates of singlet-alpha and singlet-beta product states crea-
- tes nuclear magnetisation enhancement beyond the Boltzmann le-
- vel at both the starting and the final field.
- Calculation time: minutes
- Load the spin system
- Magnet fields
- Match exchange coupling to the midpoint field
- Increase viscosity
- Get thermal equilibrium at starting field
- Set up a field ramp and time step
