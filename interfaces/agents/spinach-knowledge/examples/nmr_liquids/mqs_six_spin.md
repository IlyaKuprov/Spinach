# examples/nmr_liquids/mqs_six_spin.m

- Signature: `mqs_six_spin()`

## Purpose

Multiple-quantum (MQ) NMR experiment for a coupled system of six spins. Calculation time: minutes.

## Implementation

The example models six coupled 1H spins at 14.1 T with the specified 3J–6J couplings. It selects coherence orders +6 and -1, sets longitudinal magnetisation as the initial state and transverse magnetisation as the detection state, then simulates the refocused MQ sequence. The two-dimensional signal is squared-cosine apodised, Fourier transformed, and plotted with the 1Q and 6Q axes labelled.
