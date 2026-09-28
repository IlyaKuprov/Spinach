# examples/nmr_solids/case_studies/cp_square_vs_ramp/cp_adiabatic_vs_optimcon.m

- Signature: `cp_adiabatic_vs_optimcon()`

## Purpose

1H-15N cross-polarisation experiment in the doubly rotating frame using (a) tangent-ramped adiabatic CP; (b) numerically optimised (GRAPE method) shortcut to adiabaticity. Calculation time: minutes

## Physical / mathematical content
- Compares tangent-ramped adiabatic cross-polarisation with a GRAPE-optimised shortcut in a powder simulation, monitoring the ¹⁵N Lx signal.
- GRAPE optimises the ¹H Ly and ¹⁵N Lx control waveforms from the tangent-ramp initial guess, using L-BFGS with an SNS penalty.
- Simulates the optimised pulses at the original 500 × 2 µs timing and at half duration (500 × 1 µs).

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 1H-15N cross-polarisation experiment in the doubly rotating
- frame using (a) tangent-ramped adiabatic CP; (b) numerically
- optimised (GRAPE method) shortcut to adiabaticity.
- Calculation time: minutes
- System specification
- Interactions
- Basis set
- Spinach housekeeping
- % Tangent ramp CP simulation
- Common experiment parameters
- Simulate tangent ramped amplitude CP
- Plotting -waveform
