# examples/nmr_proteins/ct_hsqc_gb1.m

- Signature: `ct_hsqc_gb1()`

## Purpose

Constant-time HSQC experiment simulation for the GB1 protein. Simulation time: hours, faster with a Tesla A100 GPU.

## Physical / mathematical content

- Imports GB1 backbone data and removes 13C spins; the sequence observes 15N and 1H and specifies 15N decoupling in F2.

## Numerical / algorithmic content

- Simulates positive and negative FIDs with `ct_hsqc`, applies squared-cosine apodisation, Fourier-transforms and combines them into a States signal, then Fourier-transforms and plots the spectrum.

## Implementation structure

- Imports `2N9K.pdb` and `2N9K.bmrb`, sets a 14.1 T field and interaction tolerances, and builds an IK-1 `sphten-liouv` basis.
- Sets J to 90, sweeps to [3000 3000] Hz, offsets to [-7300 5100] Hz, acquisition points to [128 128], and zero filling to [512 512]; plots in ppm.
