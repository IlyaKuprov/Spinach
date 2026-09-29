# examples/nmr_liquids/pa_strychnine.m

- Signature: `pa_strychnine()`
- Source: [examples/nmr_liquids/pa_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_strychnine.m)

## What it models

A liquid-state 1H pulse-acquire spectrum for the strychnine spin system. The wrapper obtains its spin-system properties from `strychnine({'1H'})`; it does not load a measured spectrum. Its source comment describes the Redfield-superoperator treatment as an accurate line-width model and estimates a run time of seconds. Those are source comments, not results measured in this review.

## Spin system and relaxation

The field is set to 14.1 (the source labels this as the magnetic field but does not state a unit). The basis uses the spherical-tensor Liouville formalism, IK-2 approximation, scalar-coupling connectivity and proximity level 1. Redfield relaxation is enabled with zero equilibrium, retained terms set to `kite`, and a correlation time of 200e-12 s (200 ps). The greedy option is enabled and the proximity cutoff is 4.0; the wrapper gives no unit for that cutoff.

## Acquisition and processing

The initial density operator and receiver are both the 1H raising state; decoupling is empty. The wrapper calls `liquid(...,@acquire,...,'nmr')`, then applies exponential apodisation with parameter 6, Fourier transforms with 65,536 points, and plots the real spectrum. It sets offset 2800, sweep 6500, and 8,192 acquired points; those offset/sweep values and the apodisation parameter are not unit-labelled in this wrapper. The displayed axis is ppm and is inverted.

## Output and limits

The source produces an interactive 1D plot; it does not save a spectrum or report numerical peak/line-width results. The exact spin parameters supplied by the `strychnine` helper and the internals of the `acquire` callback are outside this wrapper and are not inferred here.
