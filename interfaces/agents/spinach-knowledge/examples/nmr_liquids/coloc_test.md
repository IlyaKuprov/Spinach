# examples/nmr_liquids/coloc_test.m

- Signature: `coloc_test()`

## Purpose

A simple COLOC pulse sequence example for a two-spin 1H-13C system with a long-range J-coupling. Calculation time: seconds.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A simple COLOC pulse sequence example for a two-spin
- 1H-13C system with a long-range J-coupling.
- Calculation time: seconds.
- Magnet field
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
