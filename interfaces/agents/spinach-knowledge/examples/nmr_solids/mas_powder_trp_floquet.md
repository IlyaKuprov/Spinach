# examples/nmr_solids/mas_powder_trp_floquet.m

- Signature: `mas_powder_trp_floquet()`

## Purpose

13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Isotropic chemical shifts come from the experimental data. Coordinates are from X-ray data and CSAs are estimated with DFT. Calculation time: days (hours with a Tesla card)

## Physical / mathematical content

- Models the `13C` MAS spectrum of tryptophan powder assuming `1H` decoupling, using Floquet MAS dynamics.
- Uses experimental isotropic chemical shifts, X-ray coordinates, and DFT-estimated chemical-shift anisotropies, as stated in the example description.

## Numerical / algorithmic content

- Computes Floquet FIDs for the two unit-cell molecules and sums them before exponential apodisation and Fourier transformation.

## Implementation structure

- Imports the tryptophan spin system from `trp_xray.out`, mapping C and N to `13C` and `15N`, and sets the field to 9.4 T.
- For the first molecule, applies the listed experimental isotropic shifts to the DFT Zeeman tensors; uses an `sphten-liouv`/`IK-0` basis with longitudinal `15N`, `+1` projections, and interaction level 3.
- Configures a 14 kHz rotor rate, `[1 1 1]` axis, maximum rank 11 with `leb_2ang_rank_11`, 100 kHz sweep, 2048 points, and 8192-point zero filling; computes the first molecule's FID with `floquet`.
- Reimports the spin system for a second molecule and applies its listed shifts, disables trajectory-level output, and adds its Floquet FID to the first.
- Applies exponential apodisation with parameter 6, Fourier transforms the summed FID, and plots the real spectrum.
