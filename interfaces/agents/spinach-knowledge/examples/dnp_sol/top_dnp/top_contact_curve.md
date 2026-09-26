# examples/dnp_sol/top_dnp/top_contact_curve.m

- Signature: `top_contact_curve()`

## Purpose

The transformation of -E_z into I_z during the contact time of the time-optimised pulsed DNP experiment. Further information in: Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- The transformation of -E_z into I_z during the contact time of the
- time-optimised pulsed DNP experiment. Further information in:
- Calculation time: seconds
- Q-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Basis set
- Spinach housekeeping
- Detection state
- Experiment parameters
