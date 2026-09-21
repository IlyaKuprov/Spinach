# examples/imaging/gradient_echo_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/gradient_echo_1d.m`
- Signature: `gradient_echo_1d()`
- Total lines: 64

## Purpose

A gradient echo experiment in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A gradient echo experiment in the presence
- of diffusion and flow.
- Calculation time: seconds.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Basis set
- Spinach housekeeping
- Sequence parameters
- Sample geometry

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `imaging()`, `kfigure()`, `kxlabel()`, `kylabel()`.
