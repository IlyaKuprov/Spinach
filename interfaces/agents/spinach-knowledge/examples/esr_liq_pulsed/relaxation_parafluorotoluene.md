# examples/esr_liq_pulsed/relaxation_parafluorotoluene.m

- Signature: `relaxation_parafluorotoluene()`

## Purpose

X-band pulse-acquire FFT ESR simulation of a para-fluorotoluene radical, using explicit time-domain propagation and a Redfield relaxation superoperator.

## Physical / mathematical content

- The source sets a 0.33 T field and uses secular Redfield relaxation with a correlation time of `1e-10` s.

## Numerical / algorithmic content

- The FID is acquired at 1024 points and zero-filled to 4096 points before Fourier transformation; no apodisation is applied.
- The source estimates a calculation time of hours and notes a memory requirement of at least 16 GB per CPU core.
