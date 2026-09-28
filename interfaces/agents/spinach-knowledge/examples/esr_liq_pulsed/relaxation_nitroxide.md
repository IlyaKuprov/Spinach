# examples/esr_liq_pulsed/relaxation_nitroxide.m

- Signature: `relaxation_nitroxide()`

## Purpose

Pulse-acquire FFT ESR simulation of a nitroxide radical at W-band, using explicit time-domain propagation and Redfield relaxation.

## Physical / mathematical content

- The function imports the spin system from `../standard_systems/nitroxide.log`, maps the electron and `14N` spins, and sets the field to 3.5 T.
- The basis uses the `sphten-liouv` formalism without approximation. Relaxation is Redfield with secular terms and a correlation time of `5e-11` s.

## Numerical / algorithmic content

- It acquires 512 points over a `2e8` sweep with an offset of `-2e8`, applies no apodisation, and zero-fills to 1024 points before the FFT.
- The real spectrum is plotted in the ESR setup specified by the function.
