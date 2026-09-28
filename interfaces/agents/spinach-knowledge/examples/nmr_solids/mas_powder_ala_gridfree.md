# examples/nmr_solids/mas_powder_ala_gridfree.m

- Signature: `mas_powder_ala_gridfree()`

## Purpose

Calculates the `13C` MAS spectrum of alanine powder, assuming `1H` decoupling, with a grid-free Fokker–Planck MAS formalism. The source states that magnetic parameters come from a DFT calculation and estimates seconds on a Tesla V100 GPU, much longer on CPU.

## Physical and numerical content

The spin system is constructed from the PCM-DFT alanine data in `../standard_systems/alanine.log` at 14.1 T. The basis selects the `15N` longitudinal subspace with projection +1. The calculation sets the trajectory-level option disabled and uses a 2 kHz MAS rate about `[1 1 1]`, maximum rank 17, and a 50 kHz sweep. Acquisition is for `13C` with an empty decoupling list.

## Implementation

The function calls `gridfree` with `@acquire`, applies exponential apodisation (6), zero-fills the 256-point FID to 1024 points, Fourier transforms, and plots the real spectrum. Offset is 15 kHz. A GPU-enabling line is present but commented out; the source does not enable it in the active configuration.
