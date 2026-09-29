# examples/nmr_proteins/hncoca_ubiquitin.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hncoca_ubiquitin.m)

- Signature: `hncoca_ubiquitin()`

## Purpose

A theoretical 3D HN(CO)CA simulation for human ubiquitin. The source assumes only the backbone is 13C,15N-labelled and estimates minutes of calculation time, faster with a Tesla A100 GPU. The output is simulated; the script does not load a measured 3D spectrum or establish experimental agreement.

## Protein, spin system, and acquisition

The code imports `1D3Z.pdb` and `1D3Z.bmrb` through `protein`, with molecule 1, `noshift='delete'`, and the `backbone-minimal` selection. The field is 14.1 T. Interaction/proximity cutoffs are 2.0/4.0; their units are not stated in the assignments. The basis uses `sphten-liouv`, IK-1, scalar-coupling connectivity, and interaction/proximity levels 4/1. It enables `greedy` and disables `krylov`.

The 3D sequence simulation calls `liquid(...,@hncoca,...,'nmr')` with declared spins `15N`, `13C`, and `1H`. The delay values are [2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3] s. Offsets are [-7100, 8450, 4850], sweeps [2500, 4500, 3000], acquired points [64, 64, 64], and zero-fill sizes [256, 256, 256]. Axes are displayed in ppm.

## Processing and output

Each of the four phase-cycle components receives squared-cosine apodisation. The code Fourier-transforms the four components in F3, combines conjugate components into positive and negative signals, transforms those in F2 and combines them again, then Fourier-transforms F1. The plotted 3D output is `-real(spectrum)`.
