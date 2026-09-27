# examples/nmr_liquids/hetcor_strychnine.m

- Signature: `hetcor_strychnine()`

## Purpose

Simulates a HETCOR spectrum of strychnine at natural 13C abundance. Calculation time: minutes.

## Implementation

- Loads the `1H`/`13C` strychnine spin system and sets the magnetic field to 5.9.
- Uses the `sphten-liouv` formalism with `IK-2` approximation, scalar-coupling connectivity, and proximity level 1; enables `greedy` and sets the proximity cutoff to 4.0.
- Sets `J=140`, sweep widths `[3000 10000]`, offsets `[1000 4000]`, `[256 256]` points, and `[512 512]` zero filling; decouples `1H` and uses ppm axes.
- Generates 13C isotopomers with `dilute`, simulates each with `liquid(...,@hetcor,...)` in a `parfor` loop, applies cosine apodisation in both dimensions, and sums the shifted 2D Fourier transforms before plotting the absolute spectrum.