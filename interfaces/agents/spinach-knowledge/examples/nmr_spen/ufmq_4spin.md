# examples/nmr_spen/ufmq_4spin.m

- Signature: `ufmq_4spin()`

## Purpose

Simulates a 4Q ultrafast MaxQ NMR spectrum for four coupled spins with diffusion. The source estimates hours of calculation, much faster on GPU. Authors: Maria Grazia Concilio, Ilya Kuprov, and Jean-Nicolas Dumez.

## Model and sequence

The 14.1 T system has four 1H spins with shifts 0.50, 0.35, 0.15, and 0. The listed 3J couplings are 8.0 Hz between successive spins; the listed 4J couplings are 3.0 Hz for pairs (1,3) and (2,4); the 5J coupling between spins 1 and 4 is 2.0 Hz. The example selects coherence order +4 and uses an untruncated sphten-liouv basis.

The sample length is 0.015 m with 500 points, zero flow, and diffusion coefficient 18e-10 m^2/s. Initial and detection phantoms are uniform. The imaging call uses ufmq with 120 points, 50 loops, a 6e-6 s dwell time, and a 0.041 s delay; the acquisition gradient is computed from the maximum k value. Encoding uses 500 pulse points, 40 WURST cycles, Te=0.015 s, BW=15000 Hz, Ge=0.023 T/m, and a WURST chirp.

## Processing

The example plots the k-space echo data, Fourier transforms the conventional dimension, and plots the resulting spectrum in ppm.
