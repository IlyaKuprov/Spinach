# examples/nmr_solids/mas_powder_gly_fplanck.m

- Signature: `mas_powder_gly_fplanck()`

## Purpose

Calculates the glycine powder `13C` MAS spectrum. The source header describes Fokker–Planck MAS formalism, while the implementation calls `singlerot(...)`. The header assumes `1H` decoupling; the script does not set a `parameters.decouple` value. Calculation time: seconds.

## Physical / mathematical content

- The spin system is generated from the glycine PCM-DFT log with `g2spinach`; the field is 14.1 T and the observed spin is `13C`.
- The basis uses no approximation, projection +1, and a longitudinal `15N` subspace. Interaction and proximity cutoffs are set to 5.0 and 4.0.

## Numerical / algorithmic content

- The script sets a 2000 Hz rotor rate, axis `[1 1 1]`, maximum rank 23, and grid `leb_2ang_rank_23`.
- Acquisition uses 256 points over a `5e4` sweep, zero-filled to 1024 with offset 17000; exponential apodisation parameter 6 is applied before Fourier transformation.

## Implementation structure

- Parse the glycine DFT log, create the Spinach system and basis, configure acquisition, call `singlerot(spin_system,@acquire,parameters,'nmr')`, apodise, Fourier transform, and plot.
