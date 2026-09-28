# examples/esr_liq_pulsed/relaxation_fremysalt.m

- Signature: `relaxation_fremysalt()`

## Purpose

Pulse-acquire FFT ESR simulation of Fremy salt, adapted from an EasySpin test file with acknowledgements to Stefan Stoll. The example is set to reproduce Figure 3a of [the cited paper](https://doi.org/10.1209/epl/i2004-10459-y).

## Physical / mathematical content

- The spin system contains an electron and `14N`, with anisotropic electron Zeeman and electron–nitrogen coupling tensors.
- The simulation uses explicit time propagation with a Redfield relaxation superoperator, secular terms, and a correlation time of `8e-10` s.

## Numerical / algorithmic content

- At a 0.33 T field, the function acquires an ESR FID with 512 points, applies no apodisation, and Fourier-transforms with 1024-point zero filling.
- It plots the real spectrum with first-derivative display and an inverted frequency axis.
