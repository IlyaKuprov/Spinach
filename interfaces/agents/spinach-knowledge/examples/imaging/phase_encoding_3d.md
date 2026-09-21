# examples/imaging/phase_encoding_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/phase_encoding_3d.m`
- Signature: `phase_encoding_3d()`
- Total lines: 123

## Purpose

Slice selection in 3D followed by phase-encoded imaging of the resulting slice. Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Slice selection in 3D followed by phase-encoded imaging
- of the resulting slice.
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation theory
- Disable path tracing
- This needs a GPU
- Basis set
- Spinach housekeeping
- Gat phantom from library

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `phantoms()`, `rlx_t1_t2()`, `state()`, `pulse_shape()`, `kfigure()`, `dims()`, `volplot()`, `ktitle()`, `imaging()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `apodisation()`, `fftshift()`.
