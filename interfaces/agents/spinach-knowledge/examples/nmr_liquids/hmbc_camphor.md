# examples/nmr_liquids/hmbc_camphor.m

- Signature: `hmbc_camphor()`

## Purpose

Simulates a liquid-state HMBC spectrum of camphor at natural 13C abundance using coordinates, shielding anisotropies, and J-couplings from a vacuum DFT calculation. The source notes a calculation time of seconds.

## Implementation

- Imports the camphor spin system from `../standard_systems/camphor.log` using `gparse` and `g2spinach`; sets the magnetic field to 14.1 T.
- Uses the `greedy` option, a proximity cutoff of 4.0, and an IK-2 basis with `scalar_couplings` connectivity and proximity level 1.
- Sets `J=140`, `delta_b=60e-3`, sweeps of `[40000 1500]` Hz, offsets of `[18000 900]` Hz, a `[128 128]` point grid, and `[512 512]` zero filling for `{'13C','1H'}`.
- Generates 13C isotopomers with `dilute`, simulates each with `liquid(...,@hmbc,...,'nmr')` in a `parfor` loop, applies cosine apodisation in both dimensions, and sums the shifted 2D Fourier transforms before plotting the magnitude spectrum.