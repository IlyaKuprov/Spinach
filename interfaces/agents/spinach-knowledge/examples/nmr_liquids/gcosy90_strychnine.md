# examples/nmr_liquids/gcosy90_strychnine.m

- Signature: `gcosy90_strychnine()`

## Purpose

Gradient-selected COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Gradient-selected COSY spectrum of strychnine.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform
- Form echo/anti-echo signal
