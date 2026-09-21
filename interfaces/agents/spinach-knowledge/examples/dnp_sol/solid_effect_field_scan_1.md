# examples/dnp_sol/solid_effect_field_scan_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/solid_effect_field_scan_1.m`
- Signature: `solid_effect_field_scan_1()`
- Total lines: 72

## Purpose

Magnetic field sweep DNP experiment involving a gadolinium ion, steady-state polarisation of a 15N nucleus is computed as a function of the magnetic field offset. Calculation time: minutes

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Magnetic field sweep DNP experiment involving a gadolinium ion,
- steady-state polarisation of a 15N nucleus is computed as a
- function of the magnetic field offset.
- Calculation time: minutes
- Magnetic field
- Spin system
- Electron g-tensor
- Electron ZFS tensor
- Coordinates (Angstrom)
- Basis set
- Relaxation theory
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`, `kylabel()`.
