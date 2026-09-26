# examples/dnp_sol/tppm_dnp/tppm_field_profile.m

- Signature: `tppm_field_profile()`

## Purpose

Field profile of a TPPM DNP experiment. <I_z> after a fixed contact time is calculated as a function of electron pulse amplitude and offset. Further information in: Calculation time: minutes (a large powder grid is needed)

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- Field profile of a TPPM DNP experiment. <I_z> after a fixed
- contact time is calculated as a function of electron pulse
- amplitude and offset. Further information in:
- Calculation time: minutes (a large powder grid is needed)
- Q-band magnet
- Electron and two protons
- Zeeman interactions (g-tensor for trityl, ppm guess for 1H)
- Cartesian coordinates
- Spin temperature
- Hush the output
- Basis set
- Spinach housekeeping
