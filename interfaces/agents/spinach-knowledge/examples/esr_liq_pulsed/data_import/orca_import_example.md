# examples/esr_liq_pulsed/data_import/orca_import_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/data_import/orca_import_example.m`
- Signature: `orca_import_example()`
- Total lines: 65

## Purpose

Methyl radical simulation, ORCA import. The uncommon signal intensity pattern comes from g-HFC cross-correlation.

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Methyl radical simulation, ORCA import. The uncommon signal
- intensity pattern comes from g-HFC cross-correlation.
- System properties (vacuum DFT calculation)
- Isotopes
- Zeeman interactions
- Hyperfine couplings
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `oparse()`, `gauss2mhz()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
