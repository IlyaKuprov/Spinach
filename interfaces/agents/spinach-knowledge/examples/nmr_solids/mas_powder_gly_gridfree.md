# examples/nmr_solids/mas_powder_gly_gridfree.m

- Signature: `mas_powder_gly_gridfree()`

## Purpose

Calculates the glycine powder `13C` MAS spectrum using grid-free Fokker–Planck formalism. The spin-system parameters are generated from the glycine PCM-DFT log; the script explicitly sets the field to 14.1 T. Calculation time: minutes.

## Physical / mathematical content

- `g2spinach` reads the glycine log for `13C` and `15N`; the simulation observes `13C` and uses a longitudinal `15N` subspace.
- The basis uses no approximation and projection +1. Interaction and proximity cutoffs are 5.0 and 4.0; the header assumes `1H` decoupling and the script sets `parameters.decouple={}`.

## Numerical / algorithmic content

- Grid-free acquisition uses a 2000 Hz rotor rate, axis `[1 1 1]`, and maximum rank 23.
- The FID has 256 points over a `5e4` sweep, zero-filled to 1024 with offset 17000; exponential apodisation parameter 6 is applied before Fourier transformation.

## Implementation structure

- Parse the glycine DFT log, set the field and basis options, configure the experiment, call `gridfree(spin_system,@acquire,parameters,'nmr')`, apodise, Fourier transform, and plot.
