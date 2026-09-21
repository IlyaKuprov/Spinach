# examples/nmr_solids/mas_powder_dip_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_dip_gridfree.m`
- Signature: `mas_powder_dip_gridfree()`
- Total lines: 54

## Purpose

Powder magic angle spinning spectrum of a pair of dipole-coupled proton spins using grid-free Fokker-Planck equation formalism. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Powder magic angle spinning spectrum of a pair of dipole-coupled
- proton spins using grid-free Fokker-Planck equation formalism.
- Calculation time: minutes
- System specification
- Basis set
- Spinach housekeeping
- Parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `gridfree()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
