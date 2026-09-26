# examples/esr_sol_pulsed/endor_mims_nox_powder.m

- Signature: `endor_mims_nox_powder()`

## Purpose

Mims ENDOR simulation for a nitroxide radical powder. Ideal hard pulses are assumed. Calculation time: seconds.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Mims ENDOR simulation for a nitroxide radical powder. Ideal
- hard pulses are assumed.
- Calculation time: seconds.
- Isotopes
- Magnet field
- Interactions
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Simulation
- Crude apodisation
