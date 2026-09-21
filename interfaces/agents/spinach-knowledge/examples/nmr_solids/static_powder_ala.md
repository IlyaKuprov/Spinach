# examples/nmr_solids/static_powder_ala.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_ala.m`
- Signature: `static_powder_ala()`
- Total lines: 58

## Purpose

13C NMR spectrum of static alanine powder. Protons are assumed to be decoupled. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C NMR spectrum of static alanine powder. Protons are assumed
- to be decoupled.
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

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
