# examples/nmr_proteins/ct_cosy_gb1.m

- Signature: `ct_cosy_gb1()`
- Source: [examples/nmr_proteins/ct_cosy_gb1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/ct_cosy_gb1.m)

## Purpose

Simulates a two-dimensional constant-time COSY experiment for GB1. The source estimates minutes of simulation time and notes that a Tesla A100 GPU can make it faster.

## Physical / mathematical content

- Imports molecule 1 and backbone data from `2N9K.pdb` and `2N9K.bmrb`, with `noshift='delete'`; the latter supplies protein shift information, not an experimental FID.
- Sets the field to 14.1 T, interaction and proximity cutoffs to 1.0 and 3.0, and an IK-1 `sphten-liouv` basis with scalar-coupling connectivity, `inter_level=4`, and `prox_level=1`.
- Removes `13C` and `15N` spins because the protein is treated as unlabelled. The sequence parameters contain only `1H`; this is not a 3D triple-resonance simulation.

## Numerical / algorithmic content

The caller sets offset 2400 Hz, sweeps [9000, 9000] Hz, 256 x 256 acquisition points, and 512 x 512 zero filling, with axes in ppm. It generates an FID using `liquid(...,@ct_cosy,parameters,'nmr')`, applies squared-cosine apodisation in both dimensions, performs a two-dimensional FFT and shift, and plots the magnitude spectrum. It does not load measured spectrum data.

## Implementation structure

- Imports the PDB/BMRB protein inputs, sets the field, cutoffs, basis, and greedy algorithm option.
- Builds the spin system and basis, simulates the COSY FID, apodises and Fourier-transforms it, then plots `abs(spectrum)`.
