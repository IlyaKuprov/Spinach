# examples/nmr_solids/mas_powder_gly_floquet.m

- Signature: `mas_powder_gly_floquet()`

## Purpose

Calculates a glycine powder `13C` MAS spectrum using Floquet MAS formalism. The source header says to assume `1H` decoupling; however, the script sets `parameters.decouple={}`. Calculation time: seconds.

## Physical / mathematical content

- The spin system is generated from the glycine PCM-DFT log with `g2spinach`; the script sets the field to 14.1 T and observes `13C`.
- The basis uses no approximation, projection +1, and a longitudinal `15N` subspace. Interaction and proximity cutoffs are set to 5.0 and 4.0.

## Numerical / algorithmic content

- Floquet acquisition uses a 2000 Hz rotor rate, axis `[1 1 1]`, maximum rank 23, and grid `leb_2ang_rank_23`.
- The FID has 256 points over a `5e4` sweep, zero-filled to 1024 with offset 17000; exponential apodisation parameter 6 is applied before Fourier transformation.

## Implementation structure

- Parse the glycine DFT log and generate the spin system, set field and basis options, configure the experiment, call `floquet(...)`, apodise, Fourier transform, and plot.
