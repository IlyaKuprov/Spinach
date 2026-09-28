# examples/microfluidics/plain_flow.m

- Signature: `plain_flow()`

## Purpose

Simple flow simulation with no dynamics in the spin subspace: longitudinal magnetisation is tracked as a function of time after injection into the flow field imported from COMSOL with a diffusion term also present. The tail of the pipe has drainage terms set up using a kinetics superoperator phantom.

## Physical / mathematical content

- This example imports a COMSOL mesh and velocity field and transports longitudinal magnetisation as a concentration-like quantity. It has no spin-subspace dynamics or chemical reaction; diffusion and distal-pipe drainage are included alongside flow.

## Numerical / algorithmic content

- `meshflow` evolves the Lz signal initialized in cells 140–160, using the imported mesh flow, diffusion coefficient 1e-7 m^2/s, and a distal drainage term; each trajectory frame is plotted on the mesh.

## Implementation structure

- Simple flow simulation with no dynamics in the spin subspace:
- longitudinal magnetisation is tracked as a function of time after injection into the flow field imported from COMSOL with a
- diffusion term also present. The tail of the pipe has drainage
- terms set up using a kinetics superoperator phantom.
- Import hydrodynamics information
- One proton
- Chemical shift (water)
- Basis set
- Algorithmic switches
- Spinach housekeeping
- Initial condition: Lz in a few cells
