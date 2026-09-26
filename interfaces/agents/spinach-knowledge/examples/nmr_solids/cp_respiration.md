# examples/nmr_solids/cp_respiration.m

- Signature: `cp_respiration()`

## Purpose

1H-13C RESPIRATION-CP experiment in the doubly rotating frame. Magic angle spinning simulation using Fokker-Planck formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-13C RESPIRATION-CP experiment in the doubly rotating frame.
- Magic angle spinning simulation using Fokker-Planck formalism.
- Calculation time: seconds
- Magnet field
- System specification
- Formalism and basis
- Spinach housekeeping
- Experiment parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting
