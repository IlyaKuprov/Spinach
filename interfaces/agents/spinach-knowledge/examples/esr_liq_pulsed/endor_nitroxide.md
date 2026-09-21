# examples/esr_liq_pulsed/endor_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/endor_nitroxide.m`
- Signature: `endor_nitroxide()`
- Total lines: 53

## Purpose

Mims ENDOR on a 15N-labelled nitroxide radical in liquid state. Magnetic parameters taken from a DFT calculation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Mims ENDOR on a 15N-labelled nitroxide radical in liquid
- state. Magnetic parameters taken from a DFT calculation.
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Read the spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Disable path tracing (small system)
- Sequence parameters
- Spinach housekeeping
- Simulation
- Crude apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
