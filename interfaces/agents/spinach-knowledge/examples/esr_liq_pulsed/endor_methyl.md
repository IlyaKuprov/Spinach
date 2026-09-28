# examples/esr_liq_pulsed/endor_methyl.m

- Signature: `endor_methyl()`

## Purpose

Simulates a Mims ENDOR spectrum of the liquid-state methyl radical using magnetic parameters from a DFT calculation. Calculation time: seconds.

## Physical / mathematical content

The model has three equivalent protons coupled to an electron, with an `S3`-symmetric basis for the proton spins. It detects the electron channel; the example does not specify a relaxation model.

## Numerical / algorithmic content

The sequence uses `tau=100e-9` s, 512 points over a 120 MHz sweep, and zero filling to 4096. It subtracts the mean, applies Kaiser apodisation (parameter 6), Fourier transforms, and plots the spectrum magnitude.

## Implementation structure

It constructs the full sphten-liouv basis, calls `liquid` with `@endor_mims`, and processes the FID before plotting against the nuclear-frequency axis.
