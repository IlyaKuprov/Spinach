# examples/nmr_solids/mas_powder_csa_gridfree.m

- Signature: `mas_powder_csa_gridfree()`

## Purpose

Calculates a powder MAS spectrum for a pair of anisotropically shielded proton spins using a grid-free Fokker–Planck formalism. The source estimates minutes.

## Physical and numerical content

The model declares two `1H` spins at 14.1 T, with shielding eigenvalue sets `[-2 -2 4]-5` and `[-1 -3 4]+5` and zero Euler angles. It uses a 500 Hz MAS rate about `[1 1 1]`, maximum rank 17, and a 20 kHz sweep. Acquisition is for `1H` with an empty decoupling list; no explicit powder grid is set in the source.

## Implementation

The function runs `gridfree` with `@acquire`, applies exponential apodisation (6), zero-fills the 512-point FID to 4096 points, Fourier transforms, and plots the real spectrum. The displayed axis is inverted and labelled in ppm.
