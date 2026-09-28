# examples/nmr_liquids/hmqc_strychnine.m

- Signature: `hmqc_strychnine()`

## Purpose

Simulates an HMQC spectrum of strychnine at natural 13C abundance. Calculation time: minutes.

## Implementation

- Loads the `13C`/`1H` strychnine spin system and sets the magnetic field to 5.9, with `greedy` enabled and a proximity cutoff of 4.0.
- Uses a `sphten-liouv` basis with `IK-2` approximation, scalar-coupling connectivity, and proximity level 1.
- Sets J to 140, sweeps to `[10000 3000]`, offsets to `[4000 1000]`, acquisition points to `[256 256]`, and zero filling to `[512 512]`; uses ppm axes and decouples `13C` in F2 and `1H` in F1.
- Generates `13C` isotopomers with `dilute`, simulates each with `liquid(...,@hmqc,...,'nmr')` in a `parfor` loop, applies cosine apodisation in both dimensions, and sums the shifted 2D Fourier transforms.
- Plots the magnitude spectrum with `plot_2d`.