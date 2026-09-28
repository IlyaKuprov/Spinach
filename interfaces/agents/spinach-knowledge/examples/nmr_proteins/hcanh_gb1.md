# examples/nmr_proteins/hcanh_gb1.m

- Signature: `hcanh_gb1()`

## Purpose

Simulated H(CA)NH spectrum of GB1 protein. It is assumed that only the backbone is 13C,15N-labelled. Calculation time: minutes, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Imports GB1 data with the `backbone-minimal` selection for a three-dimensional 1H, 15N, 1H experiment.

## Numerical / algorithmic content

- Simulates four FID components with `hcanh`, applies squared-cosine apodisation, Fourier-transforms F3 and F2 while combining positive and negative components, then Fourier-transforms F1 and plots the real spectrum.

## Implementation structure

- Imports `2N9K.pdb` and `2N9K.bmrb`, sets a 14.1 T field and interaction tolerances, and builds an IK-1 `sphten-liouv` basis.
- Sets sweeps to [6000 3000 6000] Hz, offsets to [4200 -7200 4200] Hz, acquisition points to [128 128 128], and zero filling to [256 256 256]; plots in ppm.
