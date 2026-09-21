# examples/nmr_solids/static_powder_nqi_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_nqi_a.m`
- Signature: `static_powder_nqi_a()`
- Total lines: 54

## Purpose

Static quadrupolar 14N powder pattern of L-valyl-L-alanine using very large numerical orientation grid, set to reproduce Figure 5 from the paper by O'Dell and Ratcliffe: Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Static quadrupolar 14N powder pattern of L-valyl-L-alanine using
- very large numerical orientation grid, set to reproduce Figure 5
- from the paper by O'Dell and Ratcliffe:
- Calculation time: minutes
- System specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
