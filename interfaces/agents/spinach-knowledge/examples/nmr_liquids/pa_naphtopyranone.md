# examples/nmr_liquids/pa_naphtopyranone.m

- Signature: `pa_naphtopyranone()`

## Purpose

Simulates the NMR spectrum of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one. The magnetic parameters are taken from [the cited source](http://dx.doi.org/10.1016/j.saa.2010.11.015). The MATLAB example comments estimate a calculation time of seconds.

## Physical / mathematical content

- Describes the liquid-state proton spin system using chemical shifts and scalar couplings from the cited magnetic parameters.
- Uses an IK-2 Liouville-space basis for the simulation.

## Numerical / algorithmic content

- Simulates a liquid-state acquisition, apodises the FID, Fourier-transforms it, and plots the spectrum.

## Implementation structure

- Defines the spin-system parameters and IK-2 basis, then runs acquisition and the apodisation/transform/plot sequence.
