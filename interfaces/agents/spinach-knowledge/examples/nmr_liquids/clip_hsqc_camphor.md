# examples/nmr_liquids/clip_hsqc_camphor.m

- Signature: `clip_hsqc_camphor()`

## Purpose

Simulates and plots a natural-abundance ¹³C CLIP-HSQC spectrum of camphor. The source estimates minutes of calculation time.

## Physical and numerical content

The molecular spin system is read from `../standard_systems/camphor.log` using `g2spinach` with min_j = 3.0 and no_xyz = 0; the source identifies coordinates, shielding anisotropies, and couplings as DFT-derived, then replaces isotropic shifts with experimental values. At 14.1 T it builds an IK-2 basis with scalar-coupling connectivity, proximity level 1, and a 4.0 proximity cutoff. It generates ¹³C isotopomers with `dilute` and simulates each using `liquid(...,@clip_hsqc,...,'nmr')` in a `parfor` loop.

The sequence uses J = 140, sweep [8000 1500], offset [4000 1000], 128 × 128 acquired points, and 512 × 512 zero filling on the ¹³C and ¹H axes (axis units: ppm). Cosine-squared apodisation is applied to the positive and negative FIDs; the code Fourier transforms the direct dimension, combines them as a States signal, transforms the indirect dimension, and plots the real spectrum.
