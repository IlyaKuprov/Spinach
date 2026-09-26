# examples/nmr_liquids/pa_strychnine.m

- Signature: `pa_strychnine()`

## Purpose

1H NMR spectrum of strychnine including an accurate model of line widths via Redfield superoperator. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H NMR spectrum of strychnine including an accurate model of
- line widths via Redfield superoperator.
- Calculation time: seconds
- Read spin system properties
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
