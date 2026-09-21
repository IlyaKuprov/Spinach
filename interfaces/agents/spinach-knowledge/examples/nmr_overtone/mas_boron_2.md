# examples/nmr_overtone/mas_boron_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_overtone/mas_boron_2.m`
- Signature: `mas_boron_2()`
- Total lines: 66

## Purpose

Overtone 10B magic angle spinning NMR spectrum. The sample is spinning in the JEOL direction. Parameters from Nghia Duong and Yusuke Nishiyama, realistic RF power and pulse width. Calculation time: hours

## Physical / mathematical content

- Overtone NMR examples. The important regime is excitation or detection of formally forbidden high-order transitions in quadrupolar nuclei, usually aided by MAS or Fokker-Planck treatments of periodic motion.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Overtone 10B magic angle spinning NMR spectrum. The sample is
- spinning in the JEOL direction. Parameters from Nghia Duong
- and Yusuke Nishiyama, realistic RF power and pulse width.
- Calculation time: hours
- System specification
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
- Magic angle
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `atan()`, `state()`, `operator()`, `singlerot()`, `kfigure()`, `plot_1d()`.
