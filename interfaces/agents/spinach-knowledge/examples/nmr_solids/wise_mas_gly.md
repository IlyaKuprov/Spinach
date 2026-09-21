# examples/nmr_solids/wise_mas_gly.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/wise_mas_gly.m`
- Signature: `wise_mas_gly()`
- Total lines: 82

## Purpose

WISE of alpha-glycine powder under MAS. Calculation time: hours, much faster on GPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- WISE of alpha-glycine powder under MAS.
- Calculation time: hours, much faster on GPU
- Spin system properties (PCM DFT calculation)
- 400 MHz spectrometer
- Isotropic alpha-glycine chemical shifts
- Basis set
- Ignore interactions below 200 Hz
- Use GPU arithmetic
- Spinach housekeeping
- Experiment setup
- Detection state
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `stack_2d()`.
