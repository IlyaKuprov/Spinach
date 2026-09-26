# examples/nmr_liquids/crazed_test.m

- Signature: `crazed_test()`

## Purpose

Long range intermolecular coherences predicted by Warren and co-workers Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Long range intermolecular coherences predicted by Warren and
- co-workers
- Calculation time: seconds
- Specify system parameters
- Use the complete basis set
- Spinach code
- Sequence parameters
- Thermal equilibrium state
- CRAZED simulation
- Apodisation
- Fourier transform
- Plotting
