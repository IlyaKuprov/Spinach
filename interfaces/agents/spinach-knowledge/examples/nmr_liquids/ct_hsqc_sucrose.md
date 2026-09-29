# examples/nmr_liquids/ct_hsqc_sucrose.m

- MATLAB implementation: [examples/nmr_liquids/ct_hsqc_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/ct_hsqc_sucrose.m)

This example computes a natural-abundance 13C/1H constant-time HSQC for sucrose, using a DFT-derived spin system with selected isotropic shifts replaced by fixed values. It is a model spectrum rather than an experimental dataset.

## Spin system and transfer

The starting parameters are parsed from the vacuum DFT log for sucrose and passed to g2spinach for 1H and 13C; the reference inputs are [31.8, 182.1], in that particle order. Import uses a 3.0 Hz scalar-coupling threshold and disables coordinate import. The example then replaces isotropic shift values for spins 1-19 and 24-30, in order, with 94.5, 73.4, 74.9, 71.5, 74.7, 62.4, 63.6, 106.0, 78.7, 76.3, 83.7, 64.7, 5.49, 3.63, 3.83, 3.54, 3.90, 3.90, 3.90, 3.75, 3.75, 4.29, 4.12, 3.96, 3.90, and 3.90 ppm, respectively. These overrides retain the calculated anisotropic shielding components while setting the isotropic shifts to the values identified as experimental in the source.

At 5.9 T, the liquid-state simulation enumerates 13C isotopomers and uses a sphten-liouv / IK-2 basis with scalar-coupling connectivity and proximity level 1; greedy selection uses a 4.0 proximity cutoff. Isotopomers are simulated in parallel with the phase-sensitive constant-time HSQC sequence and a 140 Hz working J coupling. The sequence reference links are CT-HSQC, 1992 <https://doi.org/10.1016/0022-2364(92)90144-V> and the second source cited by the sequence <https://doi.org/10.1007/BF00227470>.

## Acquisition and display

F1/F2 sweep widths are 3350/950 Hz, transmitter offsets 5000/1100 Hz, and the acquired matrix is 128 x 128 points; zero filling gives 512 points in each dimension. The 13C channel is configured for decoupling during F2 acquisition. The positive and negative FIDs each receive squared-cosine apodisation; their transformed components are combined as a States signal before the second Fourier transform. The plotted output is the real two-dimensional spectrum in ppm, using the negative-contour plotting option.

The source comments estimate seconds for the calculation. The example reports no peak assignments or experimental comparison; the hard-coded shift substitutions are spin-system inputs, not measurements from this example.
