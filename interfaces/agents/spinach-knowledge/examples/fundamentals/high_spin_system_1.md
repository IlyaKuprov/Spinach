# examples/fundamentals/high_spin_system_1.m

- Signature: `high_spin_system_1()`

## Purpose

Pulse-acquire NMR spectrum in a system with a hypothetical scalar coupling to a 235U nucleus. The spectral lines should be split accordingly.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Pulse-acquire NMR spectrum in a system with a hypothetical scalar
- coupling to a 235U nucleus. The spectral lines should be split
- accordingly.
- Magnet field
- Basis set
- Spin system
- Spinach housekeeping
- Pulse sequence parameters
- Simulation
- Apodization
- Fourier transform
- Plotting
