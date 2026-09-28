# examples/nmr_proteins/hncoca_ubiquitin.m

- Signature: `hncoca_ubiquitin()`

## Purpose

Theoretical HN(CO)CA of human ubiquitin, assuming only the backbone is 13C,15N-labelled. The source gives a calculation time of minutes and notes that a Tesla A100 GPU can make it faster.

## Setup and acquisition

The example imports `1D3Z.pdb` / `1D3Z.bmrb` with the backbone-minimal selection at 14.1 T. It uses inter/proximal cutoffs 2.0/4.0 and an IK-1 sphten-liouv basis with scalar-coupling connectivity and inter/proximal levels 4/1; greedy optimization is enabled and Krylov propagation disabled. Spins are `15N`, `13C`, and `1H`; delays [2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3] s; offsets [-7100, 8450, 4850]; sweeps [2500, 4500, 3000]; points [64, 64, 64]; zero-fill sizes [256, 256, 256]; axes are in ppm.

## Processing

The `liquid` simulation uses `@hncoca`. All four phase-cycle FIDs receive squared-cosine apodisation. Conjugate phase-cycle components are combined in the F3 and F2 transforms, then F1 is transformed; the real 3D spectrum is plotted with a negative sign.
