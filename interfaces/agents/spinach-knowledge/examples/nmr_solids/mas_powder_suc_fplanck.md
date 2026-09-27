# examples/nmr_solids/mas_powder_suc_fplanck.m

- Signature: `mas_powder_suc_fplanck()`

## Purpose

13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Chemical shielding tensors, J-couplings and coordinates are estimated with DFT. Calculation time: days

## Physical / mathematical content

- Simulates the `13C` MAS spectrum of sucrose powder assuming `1H` decoupling; chemical shielding tensors, J-couplings, and coordinates are described as DFT-derived.
- The example identifies its method as Fokker-Planck MAS; the implementation calls `singlerot` for signal acquisition.

## Numerical / algorithmic content

- Runs `singlerot` to acquire the FID, applies exponential apodisation with parameter 6, zero-fills to 1024 points, and plots the real Fourier-transformed spectrum.

## Implementation structure

- Imports the sucrose spin system from the PCM DFT log with `g2spinach`, selecting `13C`, then sets the field to 14.1 T.
- Uses the `sphten-liouv` basis with `IK-0`, `+1` projections, and interaction level 3; sets interaction and proximity cutoffs to 5.0 and 4.0.
- Configures MAS at 6000 Hz about `[1 1 1]`, maximum rank 23, a 50 kHz sweep, 256 points, 1024-point zero filling, and a 15000 Hz offset; selects `leb_2ang_rank_23` for the grid.
- Runs `singlerot` with `acquire` in NMR mode, applies exponential apodisation with parameter 6, Fourier transforms the FID, and plots the real spectrum.
