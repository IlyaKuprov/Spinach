# examples/nmr_overtone/mas_boron_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/mas_boron_1.m`
- Signature: `mas_boron_1()`
- Total lines: 54

## Purpose

Overtone Z-detection 10B magic angle spinning NMR spectrum. The sample is spinning in the JEOL direction. Parameters from Nghia Duong and Yusuke Nishiyama. The simulation focuses on the most intense of the five overtone spinning sidebands. Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Overtone Z-detection 10B magic angle spinning NMR spectrum.
- The sample is spinning in the JEOL direction. Parameters from
- Nghia Duong and Yusuke Nishiyama. The simulation focuses on
- the most intense of the five overtone spinning sidebands.
- Calculation time: hours
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `kfigure()`, `plot_1d()`.
