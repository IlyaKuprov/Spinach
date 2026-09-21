# examples/dnp_sol/cross_effect_freq_scan_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_sol/cross_effect_freq_scan_2.m`
- Signature: `cross_effect_freq_scan_2()`
- Total lines: 72

## Purpose

A simple TOTAPOL based Cross Effect DNP system. Set to repro- duce Figure 2a from Intensity differences are due to a different relaxation model and minor inconsistencies between the stated geometry and the interaction amplitudes used in the original paper. Electron rotating frame simulation using Nottingham DNP rela- xation theory detailed in Calculation time: seconds

## Physical / mathematical content

- Solid-state DNP examples. These files model microwave-driven electron-nuclear polarisation transfer mechanisms such as the solid effect, cross effect, NOVEL, XiX, TOP, BEAM, and TPPM variants. The mathematics combines driven spin dynamics, relaxation, powder/MAS averaging, and steady-state or transient propagation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A simple TOTAPOL based Cross Effect DNP system. Set to repro-
- duce Figure 2a from
- Intensity differences are due to a different relaxation model
- and minor inconsistencies between the stated geometry and the
- interaction amplitudes used in the original paper.
- Electron rotating frame simulation using Nottingham DNP rela-
- xation theory detailed in
- Calculation time: seconds
- Magnetic field
- Spin system
- Basis set
- Relaxation theory

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `crystal()`, `kfigure()`, `kxlabel()`, `kylabel()`.
