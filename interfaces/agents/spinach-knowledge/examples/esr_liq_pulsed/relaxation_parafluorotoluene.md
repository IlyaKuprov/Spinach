# examples/esr_liq_pulsed/relaxation_parafluorotoluene.m

- Signature: `relaxation_parafluorotoluene()`

## Purpose

X-band pulse-acquire FFT ESR spectrum of parafluorotoluene radical, simulated using explicit time-domain propagation including Redfield relaxation superoperator. Calculation time: minutes

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- X-band pulse-acquire FFT ESR spectrum of parafluorotoluene
- radical, simulated using explicit time-domain propagation
- including Redfield relaxation superoperator.
- Calculation time: minutes
- Ignore coordinate information (HFCs provided)
- Read the spin system (vacuum DFT calculation)
- Ignore small HFC anisotropies
- Magnet field
- Basis set
- Relaxation theory
- Spinach housekeeping
- Set the sequence parameters
