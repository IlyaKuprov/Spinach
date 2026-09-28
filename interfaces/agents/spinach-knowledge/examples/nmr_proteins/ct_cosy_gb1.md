# examples/nmr_proteins/ct_cosy_gb1.m

- Signature: `ct_cosy_gb1()`

## Purpose

Constant-time COSY experiment simulation for the GB1 protein. Simulation time: minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Imports GB1 backbone data and removes 13C and 15N spins, treating the protein as unlabelled for this proton experiment.

## Numerical / algorithmic content

- Simulates a two-dimensional proton COSY FID with `ct_cosy`, applies squared-cosine apodisation in both dimensions, then zero-fills, Fourier-transforms, shifts, and plots the magnitude spectrum.

## Implementation structure

- Imports `2N9K.pdb` and `2N9K.bmrb`, sets a 14.1 T field and interaction tolerances, and builds an IK-1 `sphten-liouv` basis.
- Uses a 2400 Hz offset, 9000 Hz sweeps, 256 × 256 acquisition points, and 512 × 512 zero filling; plots the result in ppm.
