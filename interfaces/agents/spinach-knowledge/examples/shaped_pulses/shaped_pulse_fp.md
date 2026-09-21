# examples/shaped_pulses/shaped_pulse_fp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_fp.m`
- Signature: `shaped_pulse_fp()`
- Total lines: 79

## Purpose

An off-resonance rectangular soft pulse simulated using the Fokker-Planck formalism. Note that the pulse frequency off- set accumulates as additional phase during the pulse in the same way as it would during a chirp. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- An off-resonance rectangular soft pulse simulated using the
- Fokker-Planck formalism. Note that the pulse frequency off-
- set accumulates as additional phase during the pulse in the
- same way as it would during a chirp.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Background Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `shaped_pulse_af()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
