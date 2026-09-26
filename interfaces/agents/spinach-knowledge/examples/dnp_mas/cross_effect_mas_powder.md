# examples/dnp_mas/cross_effect_mas_powder.m

- Signature: `cross_effect_mas_powder()`

## Purpose

A MAS DNP simulation performed as described in Fred Mentink- Vigier's paper (Spinach rotation conventions are different): Steady state DNP simulation for a powder. Calculation time: minutes

## Physical / mathematical content

- MAS DNP examples. These files model microwave-driven electron-nuclear polarisation transfer under magic-angle spinning, combining rotor-synchronised anisotropic interactions, relaxation, microwave irradiation, and powder/rotor averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A MAS DNP simulation performed as described in Fred Mentink-
- Vigier's paper (Spinach rotation conventions are different):
- Steady state DNP simulation for a powder.
- Calculation time: minutes
- Magnet field
- Spin specification
- Interactions
- Relaxation parameters
- Basis set
- Spinach housekeeping
- Experiment parameters
- Run the MAS DNP simulation
