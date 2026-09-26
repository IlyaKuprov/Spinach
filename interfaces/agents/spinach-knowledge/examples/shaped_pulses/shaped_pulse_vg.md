# examples/shaped_pulses/shaped_pulse_vg.m

- Signature: `shaped_pulse_vg()`

## Purpose

Veshtort-Griffin E1000B 90-degree selective pulse applied to a system of 31 proton spins with nearest-neighbor J co- uplings and linear coupling topology. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Veshtort-Griffin E1000B 90-degree selective pulse applied
- to a system of 31 proton spins with nearest-neighbor J co-
- uplings and linear coupling topology.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
