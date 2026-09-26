# examples/esr_liq_pulsed/relaxation_bisnitroxide.m

- Signature: `relaxation_bisnitroxide()`

## Purpose

X-band pulse-acquire FFT ESR spectrum of a bisnitroxide radical, using explicit time domain simulation with Redfield relaxation supeoperator. Parameters from https://doi.org/10.1039/C8CP06819D Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- X-band pulse-acquire FFT ESR spectrum of a bisnitroxide radical,
- using explicit time domain simulation with Redfield relaxation
- supeoperator. Parameters from https://doi.org/10.1039/C8CP06819D
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Couplings
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Sequence parameters
