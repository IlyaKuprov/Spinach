# examples/nmr_liquids/hmbc_cyprinol.m

- Signature: `hmbc_cyprinol()`

## Purpose

Simulates an HMBC spectrum of cyprinol at natural-abundance 13C. The source notes a calculation time of seconds.

## Implementation

- Loads the cyprinol spin system and sets the magnetic field to 11.7, with greedy mode and a proximity cutoff of 4.0.
- Uses the `sphten-liouv` basis with `IK-1` approximation, inter-level 3, proximity level 1, and scalar-coupling connectivity.
- Sets `J=150`, `delta_b=60e-3`, sweeps `[12000 2500]`, offsets `[5000 1250]`, 128 points and 512 zero-filled points in each dimension, and axes for `13C` and `1H` in ppm.
- Generates 13C isotopomers with `dilute`, simulates each using `liquid(...,@hmbc,...,'nmr')` in a parallel loop, applies cosine apodisation, and sums the shifted two-dimensional Fourier transforms before plotting the spectrum magnitude.