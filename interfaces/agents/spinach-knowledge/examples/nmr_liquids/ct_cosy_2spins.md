# examples/nmr_liquids/ct_cosy_2spins.m

- Signature: `ct_cosy_2spins()`

## Purpose

CT COSY spectrum for 2 spins. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- CT COSY spectrum for 2 spins.
- Calculation time: minutes
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting
