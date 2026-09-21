# examples/nmr_solids/rframe_nqi_dante.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/rframe_nqi_dante.m`
- Signature: `rframe_nqi_dante()`
- Total lines: 65

## Purpose

DANTE MAS spectrum of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The calculation accounts for the second-order quadrupolar shift and lineshape. Set to reproduce Figure 3d from Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- DANTE MAS spectrum of a single quadrupolar 14N nucleus using 1D
- Fokker-Planck equation and a spherical grid. The calculation
- accounts for the second-order quadrupolar shift and lineshape.
- Set to reproduce Figure 3d from
- Calculation time: minutes
- System specification
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
