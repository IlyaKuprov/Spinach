# examples/nmr_solids/mas_powder_trp_gridfree.m

- Signature: `mas_powder_trp_gridfree()`

## Purpose

13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism. Isotropic chemical shifts come from experimental data; coordinates and CSAs are estimated with DFT. The evolution generator uses a polyadic representation; see [the cited paper](https://doi.org/10.1126/sciadv.aaw8962) for further particulars. Calculation time: hours on a Tesla V100 GPU, much longer on CPU.

## Physical / mathematical content

- Models the `13C` MAS spectrum of tryptophan powder assuming `1H` decoupling; isotropic shifts are experimental, coordinates and CSAs are estimated with DFT.
- Uses grid-free Fokker-Planck MAS dynamics with a polyadic evolution generator; the example cites [the supporting paper](https://doi.org/10.1126/sciadv.aaw8962).

## Numerical / algorithmic content

- Computes and sums grid-free FIDs for the two unit-cell molecules, then applies exponential apodisation and Fourier transformation.

## Implementation structure

- Imports the tryptophan spin system from `trp_xray.out`, maps C and N to `13C` and `15N`, and sets the field to 9.4 T.
- Applies the first molecule's experimental isotropic shifts to the DFT Zeeman tensors, then uses an `sphten-liouv`/`IK-0` basis with longitudinal `15N`, `+1` projections, interaction level 3, and `greedy`/`polyadic` algorithms.
- Configures a 14 kHz rotor rate, `[1 1 1]` axis, NMR assumptions, maximum rank 11, 100 kHz sweep, 2048 points, and 8192-point zero filling; calculates the first molecule's FID with `gridfree`.
- Reimports the spin system for the second molecule, applies its experimental shifts, enables `greedy`, `polyadic`, and `gpu`, and adds its grid-free FID to the first.
- Applies exponential apodisation with parameter 6, Fourier transforms the summed FID, and plots the real spectrum.
