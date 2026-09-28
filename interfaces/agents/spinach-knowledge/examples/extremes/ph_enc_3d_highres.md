# examples/extremes/ph_enc_3d_highres.m

- Signature: `ph_enc_3d_highres()`

## Purpose

Slice selection in 3D followed by phase-encoded imaging of the resulting slice. This simulation fills up a system with eight H200 GPUs and 4 TB of RAM. Simulation time: you hope and pray this even starts; if it does, then hours.

## Physical / mathematical content

- The calculation models three-dimensional slice selection and phase-encoded imaging of the selected slice, using a single-proton spin system with T1/T2 relaxation.
- The phase-encoded acquisition is displayed in both k-space and reconstructed real-space image form.

## Numerical / algorithmic content

- The script constructs the 3D slice-selective sequence, calls `phase_enc_3d` for phase encoding, then plots the k-space data and reconstructed image.

## Implementation structure

- Slice selection in 3D followed by phase-encoded imaging
- of the resulting slice. This simulation fills up a sys-
- tem with eight H200 GPUs and 4 TB of RAM.
- Simulation time: you hope and pray this even starts;
- if it does, then hours.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation theory
- Disable path tracing
- This needs a GPU
- Basis set
