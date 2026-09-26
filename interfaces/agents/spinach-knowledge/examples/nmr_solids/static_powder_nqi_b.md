# examples/nmr_solids/static_powder_nqi_b.m

- Signature: `static_powder_nqi_b()`

## Purpose

Static powder 79Br NMR spectrum of potassium bromide. At least 3 quadrupolar tensors are necessary to reproduce the experimen- tal shape, likely due to a distribution of electrostatic envi- ronments in the powder. Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Static powder 79Br NMR spectrum of potassium bromide. At least
- 3 quadrupolar tensors are necessary to reproduce the experimen-
- tal shape, likely due to a distribution of electrostatic envi-
- ronments in the powder.
- Calculation time: seconds.
- Magnet field
- Spin system
- Chemical shift, ppm
- Quadrupolar coupling
- Basis set
- Algorithmic options
- Spinach housekeeping
