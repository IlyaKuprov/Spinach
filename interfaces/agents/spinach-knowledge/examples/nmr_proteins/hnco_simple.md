# examples/nmr_proteins/hnco_simple.m

Source: [examples/nmr_proteins/hnco_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hnco_simple.m)

- Signature: `hnco_simple()`

## Task and model

This is a synthetic 3D phase-sensitive HNCO backbone-sequence simulation, not a simulation of an imported protein or an imported measured spectrum. The four-spin model is `15N` (N), `13C` (CA), `1H` (H), and `13C` (C/CO). The field parameter is 14.1 T. Scalar shifts are [110, 55, 8.0, 180] ppm; nonzero scalar couplings are N-H 92 Hz, N-CA 11 Hz, N-C 15 Hz, and CA-C 55 Hz. The basis is `sphten-liouv` with approximation `none`.

The sequence function `hnco` specifies F1 = 15N, F2 = CO 13C, and F3 = 1H; its default receiver state is NH proton coherence. The four returned FID fields are phase/sign branches for States quadrature. The example sets sweeps [5000, 10000, 5000] Hz, offsets [-7200, 25000, 4800] Hz, and [63, 64, 65] acquired points, zero-filled to [255, 256, 257]. It sets sequence delays `tau` to [2.25, 14, 4] ms and enables the sequence's proton decoupling flag during F1; plot axes are in ppm.

## Processing and output

The script calls `liquid(spin_system,@hnco,parameters,'nmr')`, applies square-cosine apodisation to all four FID branches in all dimensions, then performs shifted FFTs along F3, F2, and F1. It adds conjugate branches when forming the absorptive F3 and F2 signals and displays `real(spectrum)` with `plot_3d`, threshold 10, bounds [0.2, 0.9, 0.2, 0.9], dimension 2, and positive contours. The source estimates seconds of calculation time; this was not benchmarked here.

The sequence source cites the phase-sensitive HNCO experiment at http://dx.doi.org/10.1016/0022-2364(90)90333-5 and the bidirectional-propagation method at http://dx.doi.org/10.1016/j.jmr.2014.04.002.
