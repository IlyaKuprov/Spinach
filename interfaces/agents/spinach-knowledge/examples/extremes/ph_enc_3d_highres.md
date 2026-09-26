# examples/extremes/ph_enc_3d_highres.m

- Signature: `ph_enc_3d_highres()`

## Purpose

Slice selection in 3D followed by phase-encoded imaging of the resulting slice. This simulation fills up a sys- tem with eight H200 GPUs and 4 TB of RAM. Simulation time: you hope and pray this even starts; if it does, then hours.

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

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
