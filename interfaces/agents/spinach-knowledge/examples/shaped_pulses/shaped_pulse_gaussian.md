# examples/shaped_pulses/shaped_pulse_gaussian.m

- Signature: `shaped_pulse_gaussian()`

## Purpose

Gaussian 90-degree pulse on a chain of 31 strongly coupled protons. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Gaussian 90-degree pulse on a chain of 31 strongly coupled protons.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
- Control and offset operators
- Initial state
