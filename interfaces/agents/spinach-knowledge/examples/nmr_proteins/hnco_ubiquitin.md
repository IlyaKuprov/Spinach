# examples/nmr_proteins/hnco_ubiquitin.m

- Signature: `hnco_ubiquitin()`

## Purpose

Theoretical HNCO of human ubiquitin, assuming only the backbone is 13C,15N-labelled. The source gives a calculation time of minutes and notes that a Tesla A100 GPU can make it faster.

## Setup and acquisition

The example imports ubiquitin from `1D3Z.pdb` and `1D3Z.bmrb` using the backbone-minimal selection and sets the field to 11.7395 T. It uses the IK-1 sphten-liouv basis with scalar-coupling connectivity (inter/proximal levels 4/1), enables greedy optimization, and disables Krylov propagation. The three acquisition dimensions are `15N`, `13C`, and `1H`; offsets are [-5900, 22050, 4100], sweeps [2000, 1600, 2900], points [64, 64, 64], and zero-fill sizes [256, 256, 256]. Delays are [2.25e-3, 14e-3, 4e-3] s; F1 decoupling is enabled and axes are in ppm.

## Processing

The `liquid` simulation uses `@hnco`. All four phase-cycle FIDs receive squared-cosine apodisation; the code Fourier-transforms F3, F2, and F1 and combines conjugate components to form the absorption-mode signal before plotting the real 3D spectrum.
