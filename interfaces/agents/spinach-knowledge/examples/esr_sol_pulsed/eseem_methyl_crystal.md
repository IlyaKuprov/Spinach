# examples/esr_sol_pulsed/eseem_methyl_crystal.m

- Signature: `eseem_methyl_crystal()`

## Purpose

Two-pulse X-band ESEEM spectrum of a methyl radical at a specific orien- tation relative to the lab frame. Magnetic parameters taken from a DFT calculation. Ideal pulses are assumed. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Two-pulse X-band ESEEM spectrum of a methyl radical at a specific orien-
- tation relative to the lab frame. Magnetic parameters taken from a DFT
- calculation. Ideal pulses are assumed.
- Calculation time: seconds
- System properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Plot the time domain signal
- Crude apodization
