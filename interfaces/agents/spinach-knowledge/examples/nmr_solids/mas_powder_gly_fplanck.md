# examples/nmr_solids/mas_powder_gly_fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_gly_fplanck.m`
- Signature: `mas_powder_gly_fplanck()`
- Total lines: 57

## Purpose

13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C MAS spectrum of glycine powder (assuming decoupling of 1H),
- computed using the Fokker-Planck MAS formalism.
- Calculation time: seconds
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
