# examples/imaging/echo_planar_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/echo_planar_2d.m`
- Signature: `echo_planar_2d()`
- Total lines: 90

## Purpose

Echo planar imaging example in 2D for a brain phantom. Simulation time: seconds, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Echo planar imaging example in 2D for a brain phantom.
- Simulation time: seconds, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rlx_t1_t2()`, `phantoms()`, `R1Ph()`, `R2Ph()`, `PDPh()`, `dims()`, `npts()`, `state()`, `imaging()`, `kfigure()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `ktitle()`.
