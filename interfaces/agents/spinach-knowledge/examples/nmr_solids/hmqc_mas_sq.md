# examples/nmr_solids/hmqc_mas_sq.m

- Signature: `hmqc_mas_sq()`

## Purpose

Powder magic angle spinning CN2D experiment (rotor-synchronized de- tection) on a 14N-1H spin pair using 1D Fokker-Planck equation and a spherical grid. The calculation accounts for the second-order qu- adrupolar shift and lineshape. Calculation time: hours on CPU, minutes with a Tesla V100 GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder magic angle spinning CN2D experiment (rotor-synchronized de-
- tection) on a 14N-1H spin pair using 1D Fokker-Planck equation and
- a spherical grid. The calculation accounts for the second-order qu-
- adrupolar shift and lineshape.
- Calculation time: hours on CPU, minutes with a Tesla V100 GPU.
- System specification
- Basis set
- Use GPU if present
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup
- Simulation
