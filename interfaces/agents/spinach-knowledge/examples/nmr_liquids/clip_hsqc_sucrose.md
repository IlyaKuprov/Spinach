# examples/nmr_liquids/clip_hsqc_sucrose.m

- Signature: `clip_hsqc_sucrose()`

## Purpose

Simulates and plots a natural-abundance ¹³C CLIP-HSQC spectrum of sucrose. The source estimates minutes of calculation time.

## Physical and numerical content

The molecular spin system comes from `../standard_systems/sucrose.log` via `g2spinach` with min_j = 3.0 and no_xyz = 0. Coordinates, shielding anisotropies, and couplings are DFT-derived in the source; selected isotropic shifts are replaced with experimental values. At 14.1 T, the code removes spins [20–23, 31–34], constructs an IK-2 scalar-coupling basis (proximity level 1; cutoff 4.0), and generates ¹³C isotopomers using `dilute`.

Each isotopomer is simulated with `liquid(...,@clip_hsqc,...,'nmr')` in parallel. Parameters are J = 140, sweep [8000 2000], offset [12000 2700], 256 × 256 acquired points, and 512 × 512 zero filling on the ¹³C and ¹H axes (axis units: ppm). Cosine-squared apodisation is applied to the positive and negative FIDs; the code Fourier transforms the direct dimension, combines them as a States signal, transforms the indirect dimension, and plots the real spectrum.
