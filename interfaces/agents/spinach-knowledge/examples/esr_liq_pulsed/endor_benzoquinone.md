# examples/esr_liq_pulsed/endor_benzoquinone.m

- Signature: `endor_benzoquinone()`

## Purpose

CW ENDOR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to reproduce Figure 2 in http://dx.doi.org/10.1002/mrc.1260280313 Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- CW ENDOR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to
- reproduce Figure 2 in http://dx.doi.org/10.1002/mrc.1260280313
- Calculation time: seconds
- Magnet field
- Isotopes and interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Crude apodization
- Fourier transform
- Plotting
