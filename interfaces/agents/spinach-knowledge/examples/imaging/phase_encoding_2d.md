# examples/imaging/phase_encoding_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/phase_encoding_2d.m`
- Signature: `phase_encoding_2d()`
- Total lines: 75

## Purpose

Simple phase-encoded 2D imaging example. Calculation time: seconds. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Implementation structure

- Simple phase-encoded 2D imaging example.
- Calculation time: seconds.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `rlx_t1_t2()`, `load()`, `state()`, `imaging()`, `get()`, `figure()`, `loc()`, `subplot()`, `mri_2d_plot()`, `ktitle()`.
