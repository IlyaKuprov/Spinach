# examples/nmr_spen/ufmq_4spin.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufmq_4spin.m)

- Signature: `ufmq_4spin()`
- Credits: Maria Grazia Concilio, Ilya Kuprov, and Jean-Nicolas Dumez.

## Experiment and spin model

This example simulates a 4Q ultrafast MaxQ NMR spectrum for four coupled protons with diffusion. At 14.1 T the 1H shifts are 0.50, 0.35, 0.15, and 0 ppm. The listed 3J couplings are 8.0 Hz for pairs (1,2), (2,3), and (3,4); the listed 4J couplings are 3.0 Hz for (1,3) and (2,4); and the 5J coupling for (1,4) is 2.0 Hz. Coherence order +4 is selected in the full sphten-liouv basis. No relaxation phantom or operator is supplied; flow is set to zero. The uniform initial phantom is longitudinal 1H magnetisation and detection is transverse 1H coherence. The calculation uses simulated data, not imported measurements.

## Spatial encoding and acquisition

The sample is 0.015 m long and represented by 500 spatial points, with derivative setting `parameters.deriv={'period',7}`. The diffusion coefficient is 18e-10 m^2/s. Acquisition uses the 1H channel, zero offset, 6.0e-6 s dwell, 120 points, and 50 loops. As in the two-spin example, the source computes maximum k from points divided by sample length, derives acquisition-gradient duration from dwell times points, and calculates the acquisition gradient from those values and the 1H spin factor.

Encoding uses 500 pulse points and 40 WURST cycles, with Te=0.015 s, bandwidth 15000 Hz, gradient Ge=0.023 T/m, and a WURST chirp. A 0.041 s delay is set. The simulation is produced by `imaging` with `@ufmq`.

## Signal display and limits

The imaginary part of the simulated k-space echo array is plotted against t2-point and k-space-point indices. The conventional dimension is Fourier transformed with a shift along dimension 2, and the magnitude is plotted in ppm for the two 1H channels specified for display. The source describes the calculation as hours, much faster on GPU; it gives no numeric runtime or measured signal. No DOI, imported measurement, or experimental validation is specified in this example.
