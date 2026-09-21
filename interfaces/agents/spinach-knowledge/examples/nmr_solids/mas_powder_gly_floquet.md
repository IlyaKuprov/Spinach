# examples/nmr_solids/mas_powder_gly_floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_gly_floquet.m`
- Signature: `mas_powder_gly_floquet()`
- Total lines: 61

## Purpose

13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C MAS spectrum of glycine powder (assuming decoupling of 1H),
- computed using the Floquet MAS formalism.
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

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `floquet()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
