# examples/nmr_spen/ufmq_6spin.m

- Signature: `ufmq_6spin()`

## Purpose

Simulates a 6Q ultrafast MaxQ NMR spectrum for six coupled spins with diffusion. The source estimates hours on an NVIDIA Tesla A100 and much longer on CPU. Authors: Maria Grazia Concilio, Ilya Kuprov, and Jean-Nicolas Dumez.

## Model and sequence

The 14.1 T system contains six 1H spins with shifts -1.0, -0.5, 0, +0.3, +0.7, and +0.9. Successive-spin 3J couplings are 8.00 Hz. The listed 4J couplings are 4.00 Hz for pairs (1,3), (1,4), (2,4), (3,5), (1,6), and (5,6); listed 5J couplings are 4.00 Hz for (2,5), (2,6), and (4,6); listed 6J couplings are 2.00 Hz for (1,5) and (3,6). The selected coherence order is +6, with an untruncated sphten-liouv basis.

The one-dimensional sample is 0.015 m long with 300 points; flow is zero and diffusion is 18e-10 m^2/s. Initial and detection phantoms are uniform. The ufmq imaging call uses 120 points, 50 loops, 6e-6 s dwell time, and a 0.041 s delay; the acquisition gradient is computed from the maximum k value. WURST encoding uses 500 pulse points, 40 cycles, Te=0.015 s, BW=15000 Hz, Ge=0.023 T/m, and chirptype wurst.

## Processing

The script displays the k-space echo data, Fourier transforms the conventional dimension, and plots the spectrum in ppm.
