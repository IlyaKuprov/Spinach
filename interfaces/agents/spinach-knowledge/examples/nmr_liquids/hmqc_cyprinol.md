# examples/nmr_liquids/hmqc_cyprinol.m

- Signature: `hmqc_cyprinol()`

## Purpose

Simulates an HMQC spectrum of cyprinol at natural 13C abundance. The source notes a calculation time of seconds.

## Implementation

- Loads the cyprinol spin system and sets the magnetic field to 11.7. Uses the greedy option, proximity cutoff 4.0, interaction cutoff 5.0, and an IK-1 spherical-tensor Liouville basis with scalar-coupling connectivity.
- Sets `J=150`, sweeps `[12000 2500]`, offsets `[5000 1250]`, 128 points and 512-point zero filling in each dimension, with spins `{'13C','1H'}` and ppm axes. Decouples 1H in F1 and 13C in F2.
- Generates 13C isotopomers with `dilute`, simulates each using `liquid(...,@hmqc,...,'nmr')` in a `parfor` loop, applies cosine apodisation, and sums their shifted two-dimensional Fourier transforms.
- Plots the absolute spectrum with `plot_2d`.