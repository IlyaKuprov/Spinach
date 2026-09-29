# examples/nmr_proteins/hncoca_simple.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hncoca_simple.m)

- Signature: `hncoca_simple()`

## Purpose

A minimal simulated 3D HNCOCA / HN(CO)CA pulse sequence, with a source-estimated calculation time of seconds. This is a four-spin model rather than a protein import or measured spectrum.

## Spin system and acquisition

The model orders its nuclei as `15N`, `13C`, `1H`, `13C`, labelled N, CA, H, and C. At 14.1 T, the scalar Zeeman shifts are [110, 55, 8, 180] ppm. The nonzero scalar couplings are N-H 92 Hz, N-CA 11 Hz, N-C 15 Hz, and CA-C 55 Hz. The basis uses `sphten-liouv` with no approximation.

The sequence call is `liquid(...,@hncoca,...,'nmr')`, with declared sequence spins `15N`, `13C`, and `1H`. Its delay values are [2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3] s. The source sets offsets to [-7200, 5600, 4800], sweeps to [5000, 8000, 5000], acquired points to [63, 64, 65], and zero-fill sizes to [255, 256, 257]. Axes are displayed in ppm.

## Processing and output

All four phase-cycle components receive squared-cosine apodisation. The code Fourier-transforms the four components in F3, combines conjugate components into positive and negative signals, transforms those in F2 and combines them again, then Fourier-transforms F1. It plots `-real(spectrum)` as a 3D spectrum.
