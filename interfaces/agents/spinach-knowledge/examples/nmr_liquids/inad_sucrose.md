# examples/nmr_liquids/inad_sucrose.m

- Signature: `inad_sucrose()`

## Purpose

INADEQUATE spectrum of sucrose. The sequence selects double-quantum coherence from coupled 13C pairs and converts it back for detection. Calculation time: minutes.

## Implementation

The example builds the sucrose spin system from a parsed vacuum-DFT log with `g2spinach`, then sets isotropic shielding values to experimental shifts. It uses an 11.7 T field and an IK-1 scalar-coupling basis. Pair-labelled 13C isotopomers are simulated with the INADEQUATE sequence (J = 50 Hz, with 1H decoupling); the accumulated FID is exponentially apodised and Fourier transformed with zero filling.
