# examples/microfluidics/plain_flow.m

- Signature: `plain_flow()`

## Purpose

Simple flow simulation with no dynamics in the spin subspace: longitudinal magnetisation is tracked as a function of time af- ter injection into the flow field imported from COMSOL with a diffusion term also present. The tail of the pipe has drainage terms set up using a kinetics superoperator phantom.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simple flow simulation with no dynamics in the spin subspace:
- longitudinal magnetisation is tracked as a function of time af-
- ter injection into the flow field imported from COMSOL with a
- diffusion term also present. The tail of the pipe has drainage
- terms set up using a kinetics superoperator phantom.
- Import hydrodynamics information
- One proton
- Chemical shift (water)
- Basis set
- Algorithmic switches
- Spinach housekeeping
- Initial condition: Lz in a few cells
