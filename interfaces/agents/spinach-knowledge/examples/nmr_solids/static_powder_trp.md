# examples/nmr_solids/static_powder_trp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_trp.m`
- Signature: `static_powder_trp()`
- Total lines: 73

## Purpose

13C NMR spectrum of tryptophan powder. Isotropic chemical shifts come from the experimental data. Coordinates and CSAs are estima- ted with DFT. Protons are assumed to be decoupled. Calculation time: hours

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 13C NMR spectrum of tryptophan powder. Isotropic chemical shifts
- come from the experimental data. Coordinates and CSAs are estima-
- ted with DFT. Protons are assumed to be decoupled.
- Calculation time: hours
- Spin system properties (DFT calculation)
- Magnet field
- Experimental chemical shifts
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
