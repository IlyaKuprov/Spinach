# examples/nmr_solids/mas_powder_nqi_floquet.m

- Signature: `mas_powder_nqi_floquet()`

## Purpose

Powder magic angle spinning spectrum of a single quadrupolar deuterium nucleus using Floquet theory. Perturbative correcti- ons to the rotationg frame transformation are not applied. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder magic angle spinning spectrum of a single quadrupolar
- deuterium nucleus using Floquet theory. Perturbative correcti-
- ons to the rotationg frame transformation are not applied.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform
- Plotting
