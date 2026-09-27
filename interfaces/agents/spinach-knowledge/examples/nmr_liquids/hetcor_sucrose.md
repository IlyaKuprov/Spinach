# examples/nmr_liquids/hetcor_sucrose.m

- Signature: `hetcor_sucrose()`

## Purpose

Simulates a HETCOR spectrum of sucrose at natural 13C abundance. Magnetic parameters come from a vacuum DFT calculation, with isotropic shielding values replaced by experimental shifts. Calculation time: minutes.

## Implementation

- Imports the sucrose spin system from `../standard_systems/sucrose.log`, sets a 5.9 T field, and uses a spherical-tensor Liouville basis with `IK-2` approximation and scalar-coupling connectivity.
- Sets `J=140`, 1H/13C sweep widths of `[1000 3350]`, offsets of `[1200 5000]`, `[256 256]` points, and `[512 512]` zero filling; decouples 1H and reports axes in ppm.
- Generates 13C isotopomers, simulates each with `hetcor` in a parallel loop, applies cosine apodisation, Fourier-transforms and sums the spectra, then plots the absolute-value 2D spectrum.