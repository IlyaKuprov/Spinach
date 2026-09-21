# examples/nmr_solids/fslghetcor_mas_gly.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/fslghetcor_mas_gly.m`
- Signature: `fslghetcor_mas_gly()`
- Total lines: 91

## Purpose

FSLG-HETCOR of alpha-glycine powder under MAS. Calculation time: hours on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- FSLG-HETCOR of alpha-glycine powder under MAS.
- Calculation time: hours on NVidia Tesla A100,
- much longer on CPU
- Spin system properties (PCM DFT calculation)
- 400 MHz spectrometer
- Isotropic alpha-glycine chemical shifts
- Basis set
- Ignore interactions below 200 Hz
- Use GPU arithmetic
- sys.enable={'gpu'};
- Spinach housekeeping
- Start with Lz of 1H, detect in quadrature on 13C

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
