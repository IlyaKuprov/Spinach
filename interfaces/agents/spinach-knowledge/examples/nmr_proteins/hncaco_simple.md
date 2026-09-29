# examples/nmr_proteins/hncaco_simple.m

Source: [examples/nmr_proteins/hncaco_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hncaco_simple.m)

- Signature: `hncaco_simple()`

## Task and model

This is a synthetic 3D HN(CA)CO backbone-sequence simulation on an explicit four-spin model, not imported measured data. The spins are `15N` (N), `13C` (CA), `1H` (H), and `13C` (C, the carbonyl label). The field parameter is 14.1 T. Scalar shifts are [110, 60, 7, 180] ppm in that order. Nonzero scalar couplings are N-H 92 Hz, N-CA 11 Hz, N-C 15 Hz, CA-C 55 Hz, CA-H 2 Hz, and H-C 4 Hz. It uses the `sphten-liouv` basis with approximation `none`.

The sequence function `hncaco` defines F1 as 15N, F2 as 13C (the carbonyl coherence for HN(CA)CO), and F3 as 1H. Its default receiver is the NH proton coherence; four returned FID fields represent phase/sign branches for States quadrature. The example sets J_NH = 92 Hz, T = 25 ms for indirect 15N evolution, and delta2 = 3 ms for coherence transfer. Sweeps are [3000, 8000, 2000] Hz, offsets [-6500, 25000, 4000] Hz, sampling is 64 points per dimension, and zero filling is [256, 256, 256]. Display axes are in ppm.

## Processing and output

The script calls `liquid(spin_system,@hncaco,parameters,'nmr')`, applies square-cosine apodisation to all four FID branches in all dimensions, and performs shifted FFTs in F3, F2, and F1. It combines conjugate pathways by subtraction for the absorptive F3 and F2 components, then plots `imag(spectrum)` with `plot_3d`, threshold 10, bounds [0.2, 0.9, 0.2, 0.9], dimension 2, and positive contours. The source estimates seconds of calculation time; this was not benchmarked here.

The sequence source cites the bidirectional-propagation method at http://dx.doi.org/10.1016/j.jmr.2014.04.002.
