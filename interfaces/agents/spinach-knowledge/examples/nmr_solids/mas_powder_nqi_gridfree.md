# examples/nmr_solids/mas_powder_nqi_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_nqi_gridfree.m`
- Signature: `mas_powder_nqi_gridfree()`
- Total lines: 56

## Purpose

Powder magic angle spinning spectrum of a single quadrupolar deuterium nucleus using grid-free Fokker-Planck MAS formalism. Second order corrections to the rotating frame transformation are not applied. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder magic angle spinning spectrum of a single quadrupolar
- deuterium nucleus using grid-free Fokker-Planck MAS formalism.
- Second order corrections to the rotating frame transformation
- are not applied.
- Calculation time: minutes
- System specification
- Basis set
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `gridfree()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
