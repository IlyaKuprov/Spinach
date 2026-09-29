# examples/nmr_liquids/clip_hsqc_sucrose.m

- MATLAB implementation: [examples/nmr_liquids/clip_hsqc_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/clip_hsqc_sucrose.m)

- Signature: `clip_hsqc_sucrose()`

## Purpose

Simulates a natural-abundance 13C CLIP-HSQC spectrum for sucrose. The input molecular geometry, shielding anisotropies and scalar couplings are described as vacuum DFT results; selected isotropic shifts are replaced with experimental values. The source estimates a calculation time of minutes.

## Spin system and sequence

The driver reads ../standard_systems/sucrose.log through g2spinach for 1H and 13C, with min_j = 3.0 and no_xyz = 0. It sets the field to 14.1 T. For spin numbers [1-19, 24-30], it substitutes isotropic shifts, in that order, with [94.5, 73.4, 74.9, 71.5, 74.7, 62.4, 63.6, 106.0, 78.7, 76.3, 83.7, 64.7, 5.49, 3.63, 3.83, 3.54, 3.90, 3.90, 3.90, 3.75, 3.75, 4.29, 4.12, 3.96, 3.90, 3.90] ppm. It then removes spins [20-23, 31-34], which the source describes as fast-exchanging or uncoupled.

The model uses a sphten-Liouville basis with IK-2 approximation, scalar-coupling connectivity and proximity level 1; greedy basis selection is enabled with proximity cutoff 4.0. The driver dilutes the spin system for 13C isotopomers and simulates each with the CLIP-HSQC pulse program in NMR mode. The program starts from proton longitudinal magnetisation, transfers coherence through J coupling using the supplied J value of 140 Hz, and observes proton coherence. Its positive and negative receiver components are apodised separately, then combined as a States signal (positive plus the conjugated negative component) before the indirect-dimension transform.

## Acquisition and display

The CLIP-HSQC pathway begins with proton longitudinal magnetisation and a proton 90-degree pulse, followed by J-coupling evolution and simultaneous 180-degree pulses on both nuclei. Subsequent proton/carbon pulses, indirect carbon evolution and coherence selection by gradient pathways produce the proton-detected signal. The driver runs this through liquid in NMR mode; it specifies no additional kinetic model. The parameter order is F1 = 13C and F2 = 1H. The sweep widths are [8000, 2000] Hz and offsets [12000, 2700] Hz, with [256, 256] acquired points and [512, 512] zero filling. Both axes are labelled in ppm. The code applies square-cosine apodisation in both dimensions, Fourier transforms the direct dimension, forms the States signal, transforms the indirect dimension, and accumulates the isotopomer spectra. It plots the real spectrum with negative contours.

The CLIP-HSQC sequence implementation cites [DOI 10.1016/j.jmr.2008.03.009](https://doi.org/10.1016/j.jmr.2008.03.009).

## Scope

This is a reduced model: the specified spins are removed, and the remaining system uses an approximate basis. Calculated coordinates, anisotropies and couplings are combined with the listed experimental isotropic shifts; the script does not provide an experimental overlay or a claim of quantitative agreement. Its plotted contours should therefore be read as the output of this stated model and processing, not as measured intensities.
