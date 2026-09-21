# examples/nmr_liquids/ct_cosy_three_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/ct_cosy_three_spin.m`
- Signature: `ct_cosy_three_spin()`
- Total lines: 55

## Purpose

CT-COSY of three spin system. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- CT-COSY of three spin system.
- Calculation time: seconds
- Spin system and interactions
- Basis set
- Algorithmic options
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
