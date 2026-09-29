# examples/nmr_proteins/hncaco_gb1.m

Source: [examples/nmr_proteins/hncaco_gb1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_proteins/hncaco_gb1.m)

- Signature: `hncaco_gb1()`

## Task and imported data

This example simulates a 3D HN(CA)CO backbone spectrum for GB1, with the source comment assuming that only the backbone is 13C,15N-labelled. It imports structure and chemical-shift inputs from `2N9K.pdb` and `2N9K.bmrb` through `protein`, selecting molecule 1, deleting spins without shifts, and using `backbone-minimal`. These are inputs to a Spinach simulation; the script does not load an experimental 3D spectrum or establish agreement with measured data. No paramagnetic centre or paramagnetic interaction is specified.

The field parameter is 14.1 T. The interaction and proximity cutoffs are 2.0 and 4.0; the example does not state their units. The basis uses `sphten-liouv`, approximation `IK-1`, connectivity `scalar_couplings`, interaction level 4, and proximity level 1. It enables `greedy` and disables `krylov`; the commented GPU option is not enabled. The sequence is HN(CA)CO: F1 is 15N, F2 is the carbonyl 13C coherence, and F3 is 1H, with the default receiver state on NH protons. The four FID fields are phase/sign branches for States quadrature.

Sequence parameters are J_NH = 92 Hz, T = 25 ms for the indirect 15N evolution delay, and delta2 = 3 ms for the coherence-transfer delay. Sweeps are [3000, 2500, 3000] Hz, offsets [-7200, 26500, 5100] Hz, and sampling is 64 points per dimension, zero-filled to 256 per dimension; the displayed axis units are ppm. The source estimates minutes of calculation, faster with a Tesla A100 GPU; that estimate is not a measured runtime here.

## Processing and output

The script invokes `liquid(spin_system,@hncaco,parameters,'nmr')`, applies square-cosine apodisation to each of the four FID branches in all dimensions, and Fourier transforms F3, F2, then F1. Conjugate branches are combined with subtraction to construct absorptive components; it plots `imag(spectrum)` with `plot_3d`, threshold 10, bounds [0.2, 0.9, 0.2, 0.9], dimension 2, and positive contours.

The sequence source cites the bidirectional-propagation method at http://dx.doi.org/10.1016/j.jmr.2014.04.002.
