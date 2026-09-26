# examples/dnp_sol/cross_effect_field_scan_1.m

- Signature: `cross_effect_field_scan_1()`

## Purpose

Magnetic field sweep cross effect DNP experiment -steady-state proton magnetisation under microwave iradiation as a function of the applied magnetic field. A powder average calculation. Calculation time: hours

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Magnetic field sweep cross effect DNP experiment -steady-state proton
- magnetisation under microwave iradiation as a function of the applied
- magnetic field. A powder average calculation.
- Calculation time: hours
- Magnetic field
- Spin system
- Electron g-tensors
- 14N quadrupolar tensor
- Coordinates (Angstrom)
- Hyperfine couplings
- Exchange coupling
- Basis set
