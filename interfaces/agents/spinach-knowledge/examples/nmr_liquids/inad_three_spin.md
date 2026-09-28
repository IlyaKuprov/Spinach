# examples/nmr_liquids/inad_three_spin.m

- Signature: `inad_three_spin()`

## Purpose

INADEQUATE spectrum of a three-spin system with J-coupling between two spins only. The sequence selects double-quantum coherence from coupled 13C pairs and converts it back for detection. Calculation time: seconds.

## Implementation

The model contains three 13C spins at 9.4 T, with chemical shifts 10, 15, and 20 ppm and a single non-zero coupling, J(1,2) = 55 Hz. Spinach simulates the one-dimensional INADEQUATE experiment with J = 55 Hz; an exponential window is applied before zero-filled Fourier transformation.
