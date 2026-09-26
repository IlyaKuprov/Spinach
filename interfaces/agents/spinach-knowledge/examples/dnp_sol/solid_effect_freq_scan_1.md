# examples/dnp_sol/solid_effect_freq_scan_1.m

- Signature: `solid_effect_freq_scan_1()`

## Purpose

A scan through the microwave frequency range in a steady state DNP experiment for a single 15N labelled urea mole- cule at a specific orientation and a specific distance from a single electron. Laboratory frame DNP simulation is carried out with state space restriction to four-spin orders and a Weizmann DNP relaxation superoperator accounting for T1 and T2 and di- polar relaxation processes. Single crystal calculatio

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A scan through the microwave frequency range in a steady
- state DNP experiment for a single 15N labelled urea mole-
- cule at a specific orientation and a specific distance
- from a single electron.
- Laboratory frame DNP simulation is carried out with state
- space restriction to four-spin orders and a Weizmann DNP
- relaxation superoperator accounting for T1 and T2 and di-
- polar relaxation processes. Single crystal calculation.
- Calculation time: minutes
- Magnetic field
- Spin system
- Basis set
