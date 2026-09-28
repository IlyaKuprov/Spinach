# examples/fundamentals/symmetry_3.m

- Signature: `symmetry_3()`

## Purpose

1H NMR spectrum of valine. Uses the fully symmetric irreducible representation of the S3(x)S3 group.

## Physical / mathematical content

The eight-proton model is specified at 11.7 T. Two groups of three equivalent protons, spins 3–5 and 6–8, are treated with S3 permutation symmetry in the spherical-tensor Liouville basis; the simulation produces a 1H valine spectrum.

## Numerical / algorithmic content

Acquires 8192 FID points with a 2500 Hz sweep and 1000 Hz offset, applies exponential apodisation with parameter 5, zero-fills to 65536 points, and Fourier-transforms the signal for plotting in ppm.

## Implementation structure

Creates the spin system and symmetric basis, configures and runs the liquid-state acquisition, then processes and plots the real spectrum.
