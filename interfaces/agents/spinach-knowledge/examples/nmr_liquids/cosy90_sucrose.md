# examples/nmr_liquids/cosy90_sucrose.m

- Signature: `cosy90_sucrose()`

## Purpose

COSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- COSY spectrum of sucrose (magnetic parameters computed with DFT).
- Calculation time: minutes
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting
