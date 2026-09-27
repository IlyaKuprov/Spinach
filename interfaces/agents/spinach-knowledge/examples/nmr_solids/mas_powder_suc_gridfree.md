# examples/nmr_solids/mas_powder_suc_gridfree.m

- Signature: `mas_powder_suc_gridfree()`

## Purpose

13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism. Chemical shielding tensors, J-couplings and coordinates are estimated with DFT. The evolution generator uses a polyadic representation; see [the cited paper](https://doi.org/10.1126/sciadv.aaw8962) for further particulars. Calculation time: hours on a Tesla V100 GPU, much longer on CPU.

## Physical / mathematical content

- Simulates the `13C` MAS spectrum of sucrose powder assuming `1H` decoupling; the source identifies shielding tensors, J-couplings, and coordinates as DFT-derived.
- Uses grid-free Fokker-Planck MAS dynamics with a polyadic representation of the evolution generator, as described in the cited paper.

## Numerical / algorithmic content

- Enables the `greedy` and `polyadic` algorithms for `gridfree`, then applies exponential apodisation with parameter 6 and Fourier transforms the FID after 1024-point zero filling.

## Implementation structure

- Imports the sucrose spin system from the PCM DFT log with `g2spinach`, selects `13C`, and sets the field to 14.1 T.
- Uses the `sphten-liouv` basis with `IK-0`, `+1` projections, and interaction level 3; sets interaction and proximity cutoffs to 5.0 and 4.0 and enables `greedy` and `polyadic`.
- Configures MAS about `[1 1 1]` at 6000 Hz, maximum rank 23, a 50 kHz sweep, 256 points, 1024-point zero filling, and a 15000 Hz offset.
- Runs `gridfree` with `acquire` in NMR mode, applies exponential apodisation with parameter 6, Fourier transforms the FID, and plots the real spectrum.
