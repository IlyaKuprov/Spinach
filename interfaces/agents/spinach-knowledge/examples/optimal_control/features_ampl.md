# examples/optimal_control/features_ampl.m

- Signature: `features_ampl()`

## Purpose

Illustrates amplitude profiling in phase-modulated pulse optimisation: the amplitude profile is prescribed in the script and the phase is optimised with LBFGS GRAPE. A DNS penalty on the second-derivative norm encourages smoothness. The spin ensemble comprises 100 equally spaced ¹³C offsets across ±160 ppm; the central 60 spins are targeted for maximum excitation, while the 20 spins on either side are unconstrained. Calculation time: minutes.

## Method

At 14.1 T, the script uses a 250-step pulse with 2 μs per step, an amplitude profile specified over those steps, and RF power levels spanning 15–20 kHz. It performs up to 1000 optimisation iterations, then simulates the free-induction decay, applies apodisation, and Fourier-transforms the signal.
