# examples/esr_liq_pulsed/relaxation_fremysalt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/relaxation_fremysalt.m`
- Signature: `relaxation_fremysalt()`
- Total lines: 72

## Purpose

Pulse-acquire FFT ESR version of the EasySpin Fremy salt test file, with acknowledgements to Stefan Stoll. The Spinach simulation is run using explicit time propagation in Liouville space with Redfield relaxation superoperator. Set to reproduce Figure 3a from Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Pulse-acquire FFT ESR version of the EasySpin Fremy salt test
- file, with acknowledgements to Stefan Stoll.
- The Spinach simulation is run using explicit time propagation
- in Liouville space with Redfield relaxation superoperator.
- Set to reproduce Figure 3a from
- Calculation time: seconds
- General layout
- Basis set
- Interactions
- Relaxation superoperator
- Spinach housekeeping
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
