# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_square_vs_ramp.m

- Signature: `cp_square_vs_ramp()`

## Purpose
Compares fixed-amplitude, linearly ramped, and tangent-ramped ¹H–¹⁵N cross-polarisation in a static-powder simulation. Calculation time: seconds.

## Physical / mathematical content
- Simulates static-powder ¹H–¹⁵N cross-polarisation with fixed-amplitude, linear-ramp, and tangent-ramp irradiation.
- Uses a 9.394 T field and a 1.05 Å ¹⁵N–¹H separation; each 1 ms contact period has 500 steps of 2 µs. The script plots the irradiation amplitudes and ¹⁵N transverse-signal expectation value.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame using (a) fixed amplitude CP; (b) linearly ramped CP;
- (c) tangent-ramped CP. Static powder simulation demonstra-
- ting the advantages of ramped cross-polarisation. For fur-
- ther information, see:
- Calculation time: seconds
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- Common experiment parameters
- Simulate fixed amplitude CP
