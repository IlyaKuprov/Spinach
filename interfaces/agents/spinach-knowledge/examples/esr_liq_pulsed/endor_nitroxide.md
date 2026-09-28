# examples/esr_liq_pulsed/endor_nitroxide.m

- Signature: `endor_nitroxide()`

## Purpose

Simulates Mims ENDOR of a 15N-labelled nitroxide radical in the liquid state using magnetic parameters from a DFT calculation. Calculation time: seconds.

## Physical / mathematical content

The spin-system parameters are imported from the supplied nitroxide Gaussian output, with the nitrogen isotope mapped as 15N. Coordinate data are ignored because hyperfine couplings are provided; the calculation detects the electron channel and uses the Mims ENDOR sequence.

## Numerical / algorithmic content

The sequence uses `tau=100e-9` s, 512 points over a 100 MHz sweep, and zero filling to 4096. It subtracts the mean, applies Kaiser apodisation (parameter 6), Fourier transforms, and plots the spectrum magnitude.

## Implementation structure

It parses the DFT data with `gparse` and `g2spinach`, builds the full sphten-liouv basis with path tracing disabled, calls `liquid` with `@endor_mims`, then processes and plots the FID on the nuclear-frequency axis.
