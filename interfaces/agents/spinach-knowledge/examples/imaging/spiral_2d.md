# examples/imaging/spiral_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/spiral_2d.m`
- Signature: `spiral_2d()`
- Total lines: 85

## Purpose

Spiral K-space imaging example in 2D. Calculation time: minutes. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Spiral K-space imaging example in 2D.
- Calculation time: minutes.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rlx_t1_t2()`, `load()`, `state()`, `imaging()`, `get()`, `figure()`, `loc()`, `subplot()`, `mri_2d_plot()`.
