# examples/nmr_liquids/pa_styrene.m

- Signature: `pa_styrene()`
- Source: [examples/nmr_liquids/pa_styrene.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_styrene.m)

## What it models

A simulated liquid-state 1H pulse-acquire spectrum of styrene, using a Zeeman-Hilbert-space model. The source presents it as a demonstration of parallel propagation and cites [the associated paper](https://doi.org/10.1063/1.3679656). Its “seconds” run-time estimate is a source comment, not a timing measured here. No experimental spectrum is loaded.

## Spin system

The wrapper explicitly defines eight 1H spins, field 7.046 (labelled magnetic induction; no unit is given in the wrapper), shifts [6.720, 6.720, 6.587, 6.587, 6.527, 6.063, 5.085, 4.555], and an 8-by-8 scalar-coupling matrix. The shift vector is supplied as chemical shifts and the plotted axis is ppm. The matrix is entered upper-triangular (zeros below the diagonal); its values include 7.7884, 7.4390, 17.6002, -0.5347 and -0.2279. The wrapper does not annotate coupling units. The basis is `zeeman-hilb` with approximation `none`; no relaxation model is configured.

## Acquisition and processing

The initial density operator and receiver are each the 1H raising state, with no decoupling. It calls `liquid(...,@acquire,...,'nmr')`, applies Gaussian apodisation with parameter 6, Fourier transforms the FID to 65,536 points and plots the real spectrum. Acquisition settings are offset 1700, sweep 1000, and 16,384 points; offset, sweep and apodisation units are not stated. The output axis is ppm and inverted. Although the example is framed around parallel propagation, this wrapper delegates acquisition to `acquire` and does not expose propagation internals.

## Output and limits

Only a plotted simulated spectrum is produced; the wrapper does not save data or compare the result against a measured spectrum. The listed spin parameters are model inputs, not reported measurements.
