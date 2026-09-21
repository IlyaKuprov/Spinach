# examples/shaped_pulses/shaped_pulse_chirp_af.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/shaped_pulses/shaped_pulse_chirp_af.m`
- Signature: `shaped_pulse_chirp_af()`
- Total lines: 87

## Purpose

Chirp inversion pulse using the Fokker-Planck formalism. Fewer points are required by the amplitude-frequency method than the "two points per period of the largest frequency" Nyquist-Shan- non condition would need for the {Cx,Cy} parameterised simula- tion of a chirped pulse. Calculation time: seconds

## Physical / mathematical content

- Shaped-pulse examples. These scripts demonstrate amplitude, phase, frequency, and gradient waveform design, including adiabatic sweeps, excitation profiles, and hardware-response considerations.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Chirp inversion pulse using the Fokker-Planck formalism. Fewer
- points are required by the amplitude-frequency method than the
- "two points per period of the largest frequency" Nyquist-Shan-
- non condition would need for the {Cx,Cy} parameterised simula-
- tion of a chirped pulse.
- Calculation time: seconds
- Magnetic field
- Isotopes
- Zeeman interactions
- Couplings
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `hamiltonian()`, `assume()`, `relaxation()`, `kinetics()`, `operator()`, `chirp_pulse()`, `shaped_pulse_af()`, `homospoil()`, `step()`, `acquire()`, `apodisation()`, `fftshift()`.
