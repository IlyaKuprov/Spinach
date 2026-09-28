# examples/nmr_solids/mas_powder_ala_fplanck.m

- Signature: `mas_powder_ala_fplanck()`

## Purpose

Calculates the `13C` MAS spectrum of alanine powder, assuming `1H` decoupling, with the Fokker–Planck MAS formalism and a spherical grid. The source estimates minutes.

## Physical and numerical content

The spin system and magnetic parameters are read from the PCM-DFT alanine calculation at `../standard_systems/alanine.log`; the field is 14.1 T. The basis selects the `15N` longitudinal subspace with projection +1. The MAS rate is 2 kHz about `[1 1 1]`, with maximum rank 17 and grid `rep_2ang_100pts_sph`. Acquisition is for `13C` with no decoupling channel specified.

## Implementation

The function builds the spin system with `g2spinach`, calls `singlerot` with `@acquire`, applies exponential apodisation (6), zero-fills the 256-point FID to 1024 points, Fourier transforms it, and plots the real spectrum. The configured sweep is 50 kHz and offset is 15 kHz.
