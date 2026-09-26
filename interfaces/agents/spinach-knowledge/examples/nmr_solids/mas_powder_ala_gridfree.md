# examples/nmr_solids/mas_powder_ala_gridfree.m

- Signature: `mas_powder_ala_gridfree()`

## Purpose

13C MAS spectrum of alanine powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism. All magnetic parameters are estimated from a DFT calculation. Calculation time: seconds on a Tesla V100 GPU, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C MAS spectrum of alanine powder (assuming decoupling of 1H),
- computed using the grid-free Fokker-Planck MAS formalism. All
- magnetic parameters are estimated from a DFT calculation.
- Calculation time: seconds on a Tesla V100 GPU,
- much longer on CPU
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Algorithmic options
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup
