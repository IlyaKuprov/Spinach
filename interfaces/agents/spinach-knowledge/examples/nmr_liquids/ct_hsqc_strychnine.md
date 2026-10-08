# examples/nmr_liquids/ct_hsqc_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/ct_hsqc_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/ct_hsqc_strychnine.m)

This example computes a natural-abundance 13C/1H constant-time HSQC for the strychnine spin system. It is a liquid-state, scalar-coupling simulation, not an imported experimental spectrum.

## Spin system and transfer

The system comes from the built-in strychnine model with 13C and 1H spins, at 5.9 T. The script dilutes the full model over 13C isotopomers and sums their calculated signals, so the result represents the natural-abundance ensemble rather than a uniformly 13C-labelled molecule. The basis uses the sphten-liouv formalism with IK-2 approximation, scalar-coupling connectivity and proximity level 1; greedy selection is enabled with a 4.0 proximity cutoff. Isotopomers are simulated in parallel with the phase-sensitive constant-time HSQC sequence using a working J coupling of 140 Hz. The sequence returns positive and negative States-quadrature components; 13C is configured for decoupling during F2 acquisition. The sequence implementation cites the original CT-HSQC paper <https://doi.org/10.1016/0022-2364(92)90144-V> and a second sequence reference <https://doi.org/10.1007/BF00227470>. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## Acquisition and display

The configured F1/F2 sweep widths are 10000/3000 Hz, transmitter offsets 4000/1000 Hz, and acquired matrix 256 x 256 points. Zero filling expands both dimensions to 512 points. Squared-cosine apodisation is applied separately to both quadrature FIDs; after the F1 transforms they are combined as a States signal, then Fourier transformed along F1. The plotted quantity is the real spectrum, with axes in ppm and the plot call requesting negative contours.

The source labels the expected calculation time as hours. No peak assignments or comparison with experimental data are supplied by this example, so the simulated contour plot should be read as a model output, not an experimental match.
