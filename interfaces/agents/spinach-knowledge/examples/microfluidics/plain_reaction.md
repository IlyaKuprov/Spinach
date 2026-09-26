# examples/microfluidics/plain_reaction.m

- Signature: `plain_reaction()`

## Purpose

Non-linear reaction kinetics in a situation when there is no hydrodynamics, diffusion, or spin dynamics. This is in- tended as a stepping stone to the more complicated cases in the same directory of the Spinach example set. Calculation time: seconds.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Non-linear reaction kinetics in a situation when there is
- no hydrodynamics, diffusion, or spin dynamics. This is in-
- tended as a stepping stone to the more complicated cases
- in the same directory of the Spinach example set.
- Calculation time: seconds.
- No spin system here
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent
- Kinetic time grid, 20 seconds
- Preallocate concentration trajectory
- Initial concentrations, mol/L
- Concentration dynamics
