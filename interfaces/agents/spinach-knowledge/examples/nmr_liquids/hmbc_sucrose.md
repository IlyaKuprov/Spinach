# examples/nmr_liquids/hmbc_sucrose.m

- Signature: `hmbc_sucrose()`

## Purpose

Simulates an HMBC spectrum of sucrose at natural 13C abundance using magnetic parameters computed with DFT. The source notes a calculation time of seconds.

## Implementation

- Imports vacuum-DFT spin-system data from `../standard_systems/sucrose.log`, maps H and C to `1H` and `13C`, and sets selected isotropic shielding values to experimental shifts.
- Uses a 5.9 T field, greedy algorithm, and a scalar-coupling-connected `sphten-liouv` basis with `IK-2` approximation.
- Sets `J=140`, `delta_b=60e-3`, sweeps of `[6000 2500]`, offsets of `[5000 900]`, a `[128 128]` grid, and `[512 512]` zero filling; axes are in ppm.
- Generates `13C` isotopomers with `dilute`, simulates each with `liquid(...,@hmbc,...,'nmr')` in a `parfor` loop, applies cosine apodisation, and sums their two-dimensional Fourier transforms before plotting the positive-magnitude spectrum.