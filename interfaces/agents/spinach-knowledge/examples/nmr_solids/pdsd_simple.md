# examples/nmr_solids/pdsd_simple.m

- Signature: `pdsd_simple()`

## Purpose

13C 2D PDSD spectrum of a simple test spin system. Calculation time: minutes, much faster on GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C 2D PDSD spectrum of a simple test spin system.
- Calculation time: minutes, much faster on GPU.
- Magnet field
- Isotopes
- Interactions (HCCH fragment)
- Basis set
- Algorithmic options
- Create the spin system structure
- Build the basis
- Experiment parameters
- Simulation
- Apodisation
