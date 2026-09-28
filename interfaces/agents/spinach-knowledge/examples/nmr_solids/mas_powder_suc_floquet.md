# examples/nmr_solids/mas_powder_suc_floquet.m

- Signature: `mas_powder_suc_floquet()`

## Purpose

13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Chemical shielding tensors, J-couplings and coordinates are estimated with DFT. Calculation time: days (hours with a Tesla card)

## Physical / mathematical content

- Simulates the `13C` MAS spectrum of sucrose powder assuming `1H` decoupling; the source identifies chemical shielding tensors, J-couplings, and coordinates as DFT-derived.
- Uses the Floquet MAS formalism to handle rotor-periodic spin dynamics.

## Numerical / algorithmic content

- Runs `floquet` for the FID, applies exponential apodisation with parameter 6, zero-fills to 1024 points for the Fourier transform, and plots the real spectrum.

## Implementation structure

- Imports the sucrose spin system from the PCM DFT log with `g2spinach`, selecting `13C`, then sets the field to 14.1 T.
- Uses the `sphten-liouv` basis with `IK-0`, `+1` projections, and interaction level 3; sets interaction and proximity cutoffs to 5.0 and 4.0.
- Configures MAS at 6000 Hz about `[1 1 1]`, maximum rank 23, a 50 kHz sweep, 256 points, 1024-point zero filling, and a 15000 Hz offset.
- Runs `floquet` with `acquire` in NMR mode, applies exponential apodisation with parameter 6, Fourier transforms the FID, and plots the real spectrum.
