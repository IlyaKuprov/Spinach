# examples/shaped_pulses/shaped_pulse_chirp_xy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_chirp_xy.m`
- Signature: `shaped_pulse_chirp_xy()`
- Total lines: 86

## Purpose

Chirped inversion pulse. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Chirped inversion pulse.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping
- Set up acquisition
- Pulse infrastructure
- Chirp waveform in amplitude-frequency coordinates
- Soft pulse

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `kinetics()`, `operator()`, `chirp_pulse()`, `shaped_pulse_xy()`, `homospoil()`, `step()`, `acquire()`, `apodisation()`, `fftshift()`.
