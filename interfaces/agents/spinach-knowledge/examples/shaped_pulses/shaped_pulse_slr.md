# examples/shaped_pulses/shaped_pulse_slr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_slr.m`
- Signature: `shaped_pulse_slr()`
- Total lines: 81

## Purpose

Shinnar-Le Roux band-selective 90-degree pulse on a chain of 31 strongly coupled protons. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Shinnar-Le Roux band-selective 90-degree pulse on a chain of
- 31 strongly coupled protons.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Hamiltonian and control operators
- Initial state
- Inversion SLR pulse waveform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `slr_pulse()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `subplot()`, `cumsum()`.
