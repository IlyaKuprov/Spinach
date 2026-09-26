# examples/nmr_solids/mas_powder_nqi_fplanck.m

- Signature: `mas_powder_nqi_fplanck()`

## Purpose

Powder magic angle spinning spectrum of a single quadrupolar deuterium nucleus using Fokker-Planck theory. Perturbative cor- rections to the rotationg frame transformation are not applied. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder magic angle spinning spectrum of a single quadrupolar
- deuterium nucleus using Fokker-Planck theory. Perturbative cor-
- rections to the rotationg frame transformation are not applied.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform
- Plotting
