# examples/esr_liq_pulsed/relaxation_parafluoronitrobenzene.m

- Signature: `relaxation_parafluoronitrobenzene()`

## Purpose

Pulse-acquire FFT ESR simulation of para-fluoronitrobenzene using explicit time-domain propagation and Redfield relaxation.

## Physical / mathematical content

- The system comprises an electron, `14N`, `19F`, and four protons, at a field of 0.33898 T.
- The relaxation superoperator is Redfield with secular terms and a 160 ps correlation time.

## Numerical / algorithmic content

- The function calls `liquid(...,@acquire,...,'esr')`, applies no apodisation, zero-fills the FID to 4096 points for the Fourier transform, and plots the real spectrum.
