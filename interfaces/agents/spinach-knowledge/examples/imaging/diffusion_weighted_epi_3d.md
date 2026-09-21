# examples/imaging/diffusion_weighted_epi_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/diffusion_weighted_epi_3d.m`
- Signature: `diffusion_weighted_epi_3d()`
- Total lines: 133

## Purpose

Three-dimensional echo planar imaging in the presence of realistic diffusion. Simulation time: hours, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Three-dimensional echo planar imaging in the presence
- of realistic diffusion.
- Simulation time: hours, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- Basis set
- Spinach housekeeping
- Pulse phase

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `pulse_shape()`, `phantoms()`, `rlx_t1_t2()`, `state()`, `kfigure()`, `dims()`, `volplot()`, `ktitle()`, `imaging()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `apodisation()`, `fftshift()`.
