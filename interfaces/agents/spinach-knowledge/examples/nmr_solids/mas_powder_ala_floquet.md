# examples/nmr_solids/mas_powder_ala_floquet.m

- Signature: `mas_powder_ala_floquet()`

## Purpose

Calculates the `13C` MAS spectrum of alanine powder, assuming `1H` decoupling, with the Floquet MAS formalism. The source estimates minutes.

## Physical and numerical content

The spin system and magnetic parameters are read from the PCM-DFT alanine calculation at `../standard_systems/alanine.log`; the field is set to 14.1 T. The calculation selects the `15N` longitudinal subspace with projection +1, and uses a MAS rate of 2 kHz about `[1 1 1]`, maximum rank 17, and the `rep_2ang_100pts_sph` grid. Acquisition is for `13C` with no decoupling channel specified.

## Implementation

The function constructs the spin system with `g2spinach`, runs `floquet` with `@acquire`, applies exponential apodisation (6), zero-fills the 256-point FID to 1024 points, Fourier transforms it, and plots the real spectrum. The configured sweep is 50 kHz and offset is 15 kHz.
