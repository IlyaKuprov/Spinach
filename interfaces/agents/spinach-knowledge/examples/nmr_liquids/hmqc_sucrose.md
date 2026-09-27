# examples/nmr_liquids/hmqc_sucrose.m

- Signature: `hmqc_sucrose()`

## Purpose

Simulates a liquid-state HMQC spectrum of sucrose at natural 13C abundance using magnetic parameters from a vacuum DFT calculation. The source notes a calculation time of seconds.

## Model and parameters

- Reads `../standard_systems/sucrose.log` for 1H and 13C, with a 3.0 Hz minimum coupling threshold, then sets selected isotropic shielding shifts to experimental values.
- Uses a 5.9 T magnetic field, a scalar-coupling-connected `sphten-liouv` basis with `IK-2` approximation, and a 4.0 proximity cutoff.
- Sets `J=140`, sweep widths `[3350 1000]`, offsets `[5000 1200]`, and a `[256 256]` acquisition grid zero-filled to `[512 512]`; axes are in ppm.

## Calculation

Generates 13C isotopomers with `dilute`, simulates each with `liquid(...,@hmqc,...)` in a `parfor` loop, applies cosine apodisation in both dimensions, and sums the shifted 2D Fourier transforms. Plots the magnitude spectrum with `plot_2d`.