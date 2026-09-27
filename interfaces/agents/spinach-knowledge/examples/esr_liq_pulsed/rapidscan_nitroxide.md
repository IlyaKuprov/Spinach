# examples/esr_liq_pulsed/rapidscan_nitroxide.m

- Signature: `rapidscan_nitroxide()`

## Purpose

Rapid scan ESR spectrum of a nitroxide radical. Calculation time: seconds

## Physical / mathematical content

- The spin system contains `14N` and an electron, with an anisotropic electron Zeeman matrix and nitrogen–electron coupling matrix at a centre field of 3.5 T.
- Relaxation uses the Redfield model with secular terms, a correlation time of `2e-11` s, and a temperature of 100 K.

## Numerical / algorithmic content

- The simulation uses the `sphten-liouv` formalism without basis approximation. It creates the spin system, sets its basis, and calls `rapidscan` with microwave power `2*pi*1e3`, a sweep from `-0.011` to `-0.003`, 500 steps, and a `1e-8` s timestep.
- The code plots the real part of the spectrum against magnetic induction.

## Implementation structure

- Rapid scan ESR spectrum of a nitroxide radical.
- Calculation time: seconds
- Centre field
- Spin system properties
- Simulation parameters
- Spinach housekeeping
- Experiment parameters
- Run the experiment
- Plot the result
