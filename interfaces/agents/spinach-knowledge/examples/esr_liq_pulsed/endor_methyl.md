# examples/esr_liq_pulsed/endor_methyl.m

- Signature: `endor_methyl()`

## Purpose

Mims ENDOR spectrum of a methyl radical in liquid state. Magnetic parameters taken from a DFT calculation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Mims ENDOR spectrum of a methyl radical in liquid state.
- Magnetic parameters taken from a DFT calculation.
- Calculation time: seconds
- Magnet field
- Spin system and interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Crude apodisation
- Fourier transform
- Plotting
