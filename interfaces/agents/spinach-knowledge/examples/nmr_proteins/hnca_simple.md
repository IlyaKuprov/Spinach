# examples/nmr_proteins/hnca_simple.m

Source: [examples/nmr_proteins/hnca_simple.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hnca_simple.m)

- Signature: `hnca_simple()`

## Task and model

This is an in-silico 3D HNCA backbone experiment on a minimal four-nucleus model, not a simulation of a supplied protein or an imported measured spectrum. The model contains `15N` (N), `13C` (CA), `1H` (H), and a second `13C` (C). The field parameter is 14.1 T. Scalar shifts are [110, 60, 8, 180] ppm in that spin order; nonzero scalar couplings are N-H 92 Hz, N-CA 11 Hz, N-C 15 Hz, CA-C 55 Hz, CA-H 2 Hz, and H-C 4 Hz.

The basis is `sphten-liouv` with approximation `none`. The sequence function `hnca` defines F1 as 15N, F2 as 13C (the CA coherence pathway in this sequence), and F3 as 1H; its default receiver state is the NH proton coherence. Its four returned fields are the phase/sign branches for States quadrature, not four different receiver nuclei. The example sets sweeps [2800, 5000, 3000] Hz, offsets [-7200, 8600, 5100] Hz, 64 points per dimension, and zero-fills each dimension to 256; axes are requested in ppm. The sequence implementation uses 92 Hz for J_NH and 11.5 Hz for J_NCA during transfer; the simple spin-system coupling entry for N-CA above is 11 Hz.

## Processing and output

The script calls `liquid(spin_system,@hnca,parameters,'nmr')`, applies square-cosine apodisation in all three dimensions, then performs shifted FFTs in F3, F2, and F1. It combines the four phase/sign branches by conjugation and addition to form the absorptive components. The plotted output is `-real(spectrum)`, using `plot_3d` with contour threshold 10, bounds [0.2, 0.9, 0.2, 0.9], dimension 2, and positive contours. The source estimates seconds of calculation time; this was not benchmarked here.

The sequence source cites the bidirectional-propagation method at http://dx.doi.org/10.1016/j.jmr.2014.04.002.
