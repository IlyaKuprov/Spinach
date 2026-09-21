# examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/pulse_acquire_benzoquinone.m`
- Signature: `pulse_acquire_benzoquinone()`
- Total lines: 75

## Purpose

Pulse-acquire FFT ESR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to reproduce Figure 1 in Simple common linewidth is used as a relaxation model. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Pulse-acquire FFT ESR on 2-methoxy-1,4-benzoquinone radical in
- liquid state. Set to reproduce Figure 1 in
- Simple common linewidth is used as a relaxation model.
- Calculation time: seconds
- Magnet induction
- Isotope list
- Zeeman interactions and couplings
- Relaxation theory
- Basis set
- Spinach housekeeping
- Experiment parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
