# examples/nmr_solids/static_powder_csa.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_csa.m`
- Signature: `static_powder_csa()`
- Total lines: 54

## Purpose

Two-spin static CSA powder pattern. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Two-spin static CSA powder pattern.
- Calculation time: seconds
- System specification
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
