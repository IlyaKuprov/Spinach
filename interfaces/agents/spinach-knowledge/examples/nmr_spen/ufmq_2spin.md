# examples/nmr_spen/ufmq_2spin.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_spen/ufmq_2spin.m)

- Signature: `ufmq_2spin()`
- Credits: Maria Grazia Concilio, Ilya Kuprov, and Jean-Nicolas Dumez.

## Experiment and spin model

This is a simulated 2Q ultrafast MaxQ NMR spectrum for two coupled protons with diffusion. At 14.1 T the 1H shifts are 0.50 and 0.15 ppm and the scalar coupling is 8.0 Hz. The source selects coherence order +2 and uses the full sphten-liouv basis. It specifies no relaxation phantom or relaxation operator, and its flow field is zero. The uniform initial phantom is longitudinal 1H magnetisation; detection is transverse 1H coherence. No measured signal is imported.

## Spatial encoding and acquisition

The sample is 0.015 m long with 500 spatial points and derivative setting `parameters.deriv={'period',7}`. The diffusion coefficient is 18e-10 m^2/s. Acquisition uses the 1H channel, zero offset, 6.0e-6 s dwell, 120 points, and 50 loops. The maximum k value is calculated as points divided by sample length; acquisition-gradient duration is dwell times points, and the gradient amplitude is calculated from that k value, duration, and the 1H spin factor.

The encoding uses 500 pulse points and 40 WURST cycles, with Te=0.015 s, bandwidth 15000 Hz, gradient Ge=0.023 T/m, and chirp type wurst. The source sets a 0.041 s delay and calls `imaging` with `@ufmq` to produce the simulated data.

## Signal display and limits

The imaginary part of the simulated k-space echo array is plotted against its t2-point and k-space-point indices. The conventional dimension is Fourier transformed with a shift along dimension 2; the magnitude is then plotted in ppm for the two 1H spins. The source estimates minutes on an NVIDIA Tesla A100 and much longer on CPU; this is not a reproduced timing. No DOI, imported measurement, or experimental validation is specified in the example.
