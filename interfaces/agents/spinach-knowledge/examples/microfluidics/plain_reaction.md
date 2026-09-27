# examples/microfluidics/plain_reaction.m

- Signature: `plain_reaction()`

## Purpose

Non-linear reaction kinetics in a situation when there is no hydrodynamics, diffusion, or spin dynamics. This is intended as a stepping stone to the more complicated cases in the same directory of the Spinach example set. Calculation time: seconds.

## Physical / mathematical content

- This standalone concentration model contains two competing second-order cycloaddition channels from cyclopentadiene and acrylonitrile to the exo and endo products. It has no spin system, hydrodynamics, or diffusion; solvent is included as an inert fifth component.

## Numerical / algorithmic content

- The nonlinear concentration-dependent reaction generator is stepped on a 20-second grid with 200 time steps using Spinach `step` and the `LG4` integrator; the plotted traces omit the solvent.

## Implementation structure

- Non-linear reaction kinetics in a situation when there is
- no hydrodynamics, diffusion, or spin dynamics. This is intended as a stepping stone to the more complicated cases
- in the same directory of the Spinach example set.
- Calculation time: seconds.
- No spin system here
- Rate constants, mol/(L*s)
- Cycloaddition reaction generator, including solvent
- Kinetic time grid, 20 seconds
- Preallocate concentration trajectory
- Initial concentrations, mol/L
- Concentration dynamics
