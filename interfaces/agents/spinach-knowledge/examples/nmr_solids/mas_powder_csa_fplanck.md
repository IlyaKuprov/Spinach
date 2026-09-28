# examples/nmr_solids/mas_powder_csa_fplanck.m

- Signature: `mas_powder_csa_fplanck()`

## Purpose

Calculates the powder MAS spectrum of a pair of anisotropically shielded protons using a Fokker–Planck-based formalism. The source estimates seconds.

## Physical and numerical content

The system has two `1H` spins at 14.1 T, with shielding eigenvalue sets `[-2 -2 4]-5` and `[-1 -3 4]+5` and zero Euler angles. The MAS rate is 500 Hz about `[1 1 1]`, with maximum rank 17 and grid `leb_2ang_rank_17`. Acquisition is for `1H`.

## Implementation

The function calls `singlerot` with `@acquire`, applies exponential apodisation (6), zero-fills the 512-point FID to 4096 points, Fourier transforms, and plots the real spectrum. The sweep is 20 kHz and the axis units are ppm.
