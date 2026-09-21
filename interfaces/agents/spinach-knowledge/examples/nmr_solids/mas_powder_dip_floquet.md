# examples/nmr_solids/mas_powder_dip_floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_dip_floquet.m`
- Signature: `mas_powder_dip_floquet()`
- Total lines: 56

## Purpose

Spinning powder pulse-acquire experiment on a two-spin system with a dipolar coupling using Floquet theory. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Spinning powder pulse-acquire experiment on a two-spin system
- with a dipolar coupling using Floquet theory.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Pulse-acquire setup
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `floquet()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
