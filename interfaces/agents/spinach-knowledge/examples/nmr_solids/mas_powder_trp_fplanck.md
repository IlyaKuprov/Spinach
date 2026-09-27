# examples/nmr_solids/mas_powder_trp_fplanck.m

- Signature: `mas_powder_trp_fplanck()`

## Purpose

13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Isotropic chemical shifts come from the experimental data. Coordinates are from X-ray data and CSAs are estimated with DFT. Calculation time: days, hours with a Tesla A100 GPU.

## Physical / mathematical content

- Models the `13C` MAS spectrum of tryptophan powder assuming `1H` decoupling; isotropic shifts are experimental, coordinates are from X-ray data, and CSAs are estimated with DFT.
- The source describes the calculation as Fokker-Planck MAS and uses `singlerot` to acquire each molecule's signal.

## Numerical / algorithmic content

- Acquires and sums the `singlerot` FIDs for the two unit-cell molecules, then applies exponential apodisation and a Fourier transform.

## Implementation structure

- Imports the tryptophan spin system from `trp_xray.out`, mapping C and N to `13C` and `15N`, and sets the field to 9.4 T.
- For each of two molecules in the unit cell, applies the corresponding experimental isotropic shifts to the DFT Zeeman tensors. The first uses an `sphten-liouv`/`IK-0` basis with longitudinal `15N`, `+1` projections, and interaction level 3; the second enables GPU execution.
- Configures a 14 kHz rotor rate, `[1 1 1]` axis, maximum rank 11 with `leb_2ang_rank_11`, 100 kHz sweep, 2048 points, and 8192-point zero filling.
- Runs `singlerot` with `acquire` in NMR mode for each molecule and sums the FIDs.
- Applies exponential apodisation with parameter 6, Fourier transforms the summed FID, and plots the real spectrum.
