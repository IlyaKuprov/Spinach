# examples/nmr_liquids/inv_rec_strychnine.m

- Signature: `inv_rec_strychnine()`
- Source: [`examples/nmr_liquids/inv_rec_strychnine.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/inv_rec_strychnine.m)

## Purpose

A simulated homonuclear 1H inversion-recovery experiment for strychnine. The source identifies the example as 250 MHz and gives an estimated calculation time of minutes.

## Spin system and relaxation

The model comes from `strychnine({'1H'})`; the magnetic-field parameter is 5.9. The calculation enables greedy parallelisation and uses Redfield relaxation, Di Bari equilibrium, `rlx_keep='kite'`, a correlation-time entry of `200e-12`, and a temperature parameter of 298 (the source does not annotate units for these parameter values). The basis is `sphten-liouv` with the IK-2 scalar-coupling approximation, scalar-coupling connectivity, proximity level 1, and a proximity cutoff of 5.0.

## Sequence and spectrum

The liquid-NMR simulation calls `liquid(spin_system,@inv_rec,parameters,'nmr')` for the 1H channel. It sets offset 1250, sweep 2500, 4096 points, a maximum-delay parameter of 1.0, ten delays, ppm axis units, and axis inversion. The source does not annotate units for offset, sweep, or maximum delay, nor list the individual delay values. The resulting FIDs receive exponential apodisation with parameter 6; an FFT along the acquisition dimension is shifted and the real spectrum is plotted.

The source provides no fitted T1 value, measured spectrum, or DOI. Its numerical output is a simulation, not an experimental rate determination.
