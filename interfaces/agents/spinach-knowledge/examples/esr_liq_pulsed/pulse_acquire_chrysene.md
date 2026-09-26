# examples/esr_liq_pulsed/pulse_acquire_chrysene.m

- Signature: `pulse_acquire_chrysene()`

## Purpose

W-band pulse-acquire FFT ESR spectrum of a chrysene cation radical in a non-viscous liquid. Simple common line width is used as a relaxation model. Symmetry treatment is performed using the full S2xS2xS2xS2xS2xS2 group direct product. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- W-band pulse-acquire FFT ESR spectrum of a chrysene cation radical in
- a non-viscous liquid. Simple common line width is used as a relaxation
- model. Symmetry treatment is performed using the full S2xS2xS2xS2xS2xS2
- group direct product.
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Read the spin system (vacuum DFT calculation)
- Magnet induction
- Relaxation theory
- Basis set
- Symmetry
- Spinach housekeeping
