# examples/nmr_solids/mas_powder_suc_floquet.m

- Signature: `mas_powder_suc_floquet()`

## Purpose

13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Chemical shielding tensors, J-couplings and coordinates are estimated with DFT. Calculation time: days (hours with a Tesla card)

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C MAS spectrum of sucrose powder (assuming decoupling of 1H),
- computed using the Floquet MAS formalism. Chemical shielding
- tensors, J-couplings and coordinates are estimated with DFT.
- Calculation time: days (hours with a Tesla card)
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- % Simulation
- Apodisation
