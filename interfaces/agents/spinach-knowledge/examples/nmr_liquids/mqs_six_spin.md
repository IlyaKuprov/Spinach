# examples/nmr_liquids/mqs_six_spin.m

- Signature: `mqs_six_spin()`

## Purpose

Multiple-quantum (MQ) NMR experiment for a coupled system of six spins. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Multiple-quantum (MQ) NMR experiment for a coupled system
- of six spins.
- Calculation time: minutes
- Magnetic field
- Chemical shifts
- 3J couplings
- 4J couplings
- 5J couplings
- 6J couplings
- Coherence to select
- Basis set
- Algorithmic options
