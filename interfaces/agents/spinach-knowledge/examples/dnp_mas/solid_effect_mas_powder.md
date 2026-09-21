# examples/dnp_mas/solid_effect_mas_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_mas/solid_effect_mas_powder.m`
- Signature: `solid_effect_mas_powder()`
- Total lines: 65

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
- Stack generation parameters
- Run the MAS DNP simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `masdnp()`, `num2str()`.
