# examples/nmr_spen/psycosy_acrolein.m

- Signature: `psycosy_acrolein()`

## Purpose

PSYCOSY of Acrolein. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- The source specifies a five-proton acrolein spin system at 14.1 T, including chemical shifts and scalar couplings.
- The FID is simulated with the `@psycosy` sequence through the imaging function.

## Numerical / algorithmic content

- The simulated FID is apodised with square-sine windows in both dimensions, then transformed with a 2D FFT. The plotted result is the magnitude spectrum.

## Implementation structure

- Defines the spin system, basis, and sequence parameters; runs imaging with `@psycosy`; applies apodisation and plots the 2D Fourier spectrum.
