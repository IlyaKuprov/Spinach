# examples/nmr_solids/mas_powder_csa_floquet.m

- Signature: `mas_powder_csa_floquet()`

## Purpose

The source describes a powder MAS spectrum of a single anisotropically shielded proton, using a Floquet-based formalism. It estimates seconds.

## Physical and numerical content

The implementation sets a 14.1 T field and declares two `1H` spins, with separate shielding eigenvalue sets `[-2 -2 4]-5` and `[-1 -3 4]+5` and zero Euler angles. Thus, the header's “single” proton description does not match the two-spin system declaration. The MAS rate is 500 Hz about `[1 1 1]`; the Floquet grid is `leb_2ang_rank_17`, with maximum rank 17.

## Implementation

The function runs `floquet` with `@acquire`, applies exponential apodisation (6), zero-fills the 512-point FID to 4096 points, Fourier transforms, and plots the real spectrum. The sweep is 20 kHz, with zero offset and inverted ppm axis.
