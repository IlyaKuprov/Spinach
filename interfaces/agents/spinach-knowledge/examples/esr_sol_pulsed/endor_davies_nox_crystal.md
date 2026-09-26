# examples/esr_sol_pulsed/endor_davies_nox_crystal.m

- Signature: `endor_davies_nox_crystal()`

## Purpose

Davies ENDOR simulation for a nitroxide radical at a single orientation. Soft pulses are simulated using Fokker-Planck formalism. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Davies ENDOR simulation for a nitroxide radical at a single
- orientation. Soft pulses are simulated using Fokker-Planck
- formalism.
- Calculation time: minutes
- Isotopes
- Magnet field
- Interactions
- Basis set
- Relaxation theory
- Spinach housekeeping
- % Stage 1: pulse-acquire ESR spectrum
- Sequence parameters
