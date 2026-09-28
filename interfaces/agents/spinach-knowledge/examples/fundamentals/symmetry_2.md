# examples/fundamentals/symmetry_2.m

- Signature: `symmetry_2()`

## Purpose

Pulse-acquire NMR spectrum of a highly symmetric spin system provided by Andres Castillo. Uses the fully symmetric irreducible representation of the S3(x)S3(x)S3 permutation symmetry group.

## Physical / mathematical content

The source specifies 13 protons at 9.4 T and three S3 groups acting on spins 1–3, 4–6, and 7–9. The spherical-tensor Liouville basis uses the IK-2 approximation, scalar-coupling connectivity, proximity level 1, and projection +1.

## Numerical / algorithmic content

Simulates a liquid-state pulse-acquire FID with 2048 points, a 2000 Hz sweep and 800 Hz offset; applies exponential apodisation with parameter 6, zero-fills to 8196 points, and Fourier-transforms and plots the real spectrum in ppm.

## Implementation structure

Builds the symmetric basis for the supplied spin system, configures 1H acquisition, runs the liquid simulation, and processes the resulting FID.
