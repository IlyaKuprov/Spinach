# examples/nmr_proteins/hncoca_simple.m

- Signature: `hncoca_simple()`

## Purpose

A minimal HNCOCA pulse-sequence simulation; the source estimates a calculation time of seconds.

## Spin system and acquisition

The four-spin model is ordered `15N`, `13C`, `1H`, `13C` and labelled N, CA, H, and C. The field is 14.1 T; the listed scalar Zeeman shifts are [110, 55, 8, 180], and the nonzero couplings are N–H 92, N–CA 11, N–C 15, and CA–C 55 (source values). The basis uses the sphten-liouv formalism with no approximation. Sequence delays are [2.25e-3, 2.75e-3, 8.00e-3, 7.00e-3] s; spins are `15N`, `13C`, `1H`, offsets [-7200, 5600, 4800], sweeps [5000, 8000, 5000], points [63, 64, 65], and zero-fill sizes [255, 256, 257]. Axes are in ppm.

## Processing

The `liquid` simulation uses `@hncoca`. The four phase-cycle components receive squared-cosine apodisation; conjugate components are combined during the F3 and F2 transforms, followed by the F1 transform. The plotted 3D spectrum is `-real(spectrum)`.
