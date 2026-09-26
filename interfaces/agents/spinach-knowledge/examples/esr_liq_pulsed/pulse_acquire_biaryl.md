# examples/esr_liq_pulsed/pulse_acquire_biaryl.m

- Signature: `pulse_acquire_biaryl()`

## Purpose

A time-domain pulse-acquire version of the EasySpin biaryl test file, with acknowledgements to Stefan Stoll. The Spinach simulation is run using explicit time propagation in Liou- ville space. Full symmetry treatment using the S2xS2xS2xS2xS2xS2 group direct product is performed. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A time-domain pulse-acquire version of the EasySpin biaryl test file,
- with acknowledgements to Stefan Stoll.
- The Spinach simulation is run using explicit time propagation in Liou-
- ville space. Full symmetry treatment using the S2xS2xS2xS2xS2xS2 group
- direct product is performed.
- Calculation time: seconds
- Magnet induction
- Isotope list
- Basis set
- Zeeman interactions
- Spin-spin couplings
- Relaxation theory
