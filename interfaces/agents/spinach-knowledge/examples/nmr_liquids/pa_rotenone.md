# examples/nmr_liquids/pa_rotenone.m

- Signature: `pa_rotenone()`

## Purpose

Simulates the ¹H NMR spectrum of rotenone with a T1/T2 relaxation model. Magnetic parameters are taken from [the cited source](http://dx.doi.org/10.1002/jhet.5570250160). The MATLAB example comments estimate a calculation time of seconds.

## Physical / mathematical content

- Models a liquid-state rotenone proton spin system with T1/T2 relaxation.

## Numerical / algorithmic content

- Runs liquid-state acquisition, applies Gaussian apodisation to the FID, Fourier-transforms it, and plots the real spectrum.

## Implementation structure

- Defines the 22-proton system and relaxation parameters, simulates the FID, and performs apodisation, Fourier transformation, and plotting.
