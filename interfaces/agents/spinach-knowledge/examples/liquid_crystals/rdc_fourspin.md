# examples/liquid_crystals/rdc_fourspin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/liquid_crystals/rdc_fourspin.m`
- Signature: `rdc_fourspin()`
- Total lines: 70

## Purpose

CLIP-HSQC spectrum of a four-spin system in a liquid crystal with a user-specified order matrix. Calculation times: seconds.

## Physical / mathematical content

- Liquid-crystal examples. These scripts exploit partial ordering and Saupe-tensor physics, so anisotropic couplings survive orientational averaging and generate residual dipolar couplings or anisotropic transfer behaviour.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- CLIP-HSQC spectrum of a four-spin system in a liquid crystal
- with a user-specified order matrix.
- Calculation times: seconds.
- Magnet field
- Spin system and interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
