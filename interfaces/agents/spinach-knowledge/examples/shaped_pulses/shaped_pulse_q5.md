# examples/shaped_pulses/shaped_pulse_q5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_q5.m`
- Signature: `shaped_pulse_q5()`
- Total lines: 91

## Purpose

90-degree Q5 pulse on a chain of 31 strongly coupled protons. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 90-degree Q5 pulse on a chain of 31 strongly coupled protons.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Assumptions
- Hamiltonian superoperator
- Control and offset operators
- Initial state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `read_wave()`, `polar2cartesian()`, `shaped_pulse_xy()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
