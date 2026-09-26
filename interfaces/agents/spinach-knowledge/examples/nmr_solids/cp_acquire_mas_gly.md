# examples/nmr_solids/cp_acquire_mas_gly.m

- Signature: `cp_acquire_mas_gly()`

## Purpose

1H-13C cross-polarisation followed by acquisition under magic angle spinning in alpha-glycine powder. Reduced Liouville spa- ce is used: up to, and including, three-spin correlations. Calculation time: minutes on Tesla A100, much longer on CPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H-13C cross-polarisation followed by acquisition under magic
- angle spinning in alpha-glycine powder. Reduced Liouville spa-
- ce is used: up to, and including, three-spin correlations.
- Calculation time: minutes on Tesla A100, much longer on CPU.
- Spin system properties (PCM DFT calculation)
- 400 MHz spectrometer
- Isotropic alpha-glycine chemical shifts
- Spin temperature
- Basis set
- Algorithmic options
- Neglect interactions below 200 Hz
- Spinach housekeeping
