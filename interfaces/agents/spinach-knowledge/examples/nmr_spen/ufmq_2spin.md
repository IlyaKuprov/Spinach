# examples/nmr_spen/ufmq_2spin.m

- Signature: `ufmq_2spin()`

## Purpose

Simulates a 2Q ultrafast MaxQ NMR spectrum for two coupled spins with diffusion. The source estimates minutes on an NVIDIA Tesla A100 and much longer on CPU. Authors: Maria Grazia Concilio, Ilya Kuprov, and Jean-Nicolas Dumez.

## Model and sequence

The 14.1 T system comprises two 1H spins with shifts 0.50 and 0.15 and an 8.0 Hz scalar coupling. It selects coherence order +2 and uses an untruncated sphten-liouv basis. The one-dimensional sample length is 0.015 m with 500 points; flow is zero and diffusion is 18e-10 m^2/s. Initial and detection state phantoms are uniform, with Lz and L+ 1H states.

The imaging call uses ufmq. Acquisition uses 120 points, 50 loops, 6e-6 s dwell time, zero offset, and a delay of 0.041 s; the acquisition gradient is calculated from the maximum k value. WURST encoding uses 500 pulse points, 40 WURST cycles, Te=0.015 s, BW=15000 Hz, Ge=0.023 T/m, and the wurst chirp type.

## Processing

The example plots the k-space echo data, Fourier transforms the conventional dimension, and plots the resulting spectrum in ppm.
