# examples/optimal_control/case_studies/Tosner_JMR_2009/bb_inversion_pulse.m

- Signature: `bb_inversion_pulse()`

## Purpose

Broadband inversion pulse design for liquid-state NMR. Reproduces the second example from the cited paper ([Tosner et al., JMR (2009)](http://dx.doi.org/10.1016/j.jmr.2008.11.020)) using Spinach. The source models a single 1H spin at 14.1 T in the rotating frame, with transmitter offsets. It designs a 600 µs pulse in 600 one-microsecond slices to invert Iz to −Iz over a ±50 kHz design range using Cartesian Lx and Ly controls.

## Implementation

- Optimizes over 101 design offsets from −50 to +50 kHz using GRAPE with L-BFGS, a 10 kHz power level, and a 200-iteration maximum.
- Simulates inversion efficiency over 201 offsets from −100 to +100 kHz and plots the resulting offset profile.
