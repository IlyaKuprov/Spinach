# examples/nmr_liquids/cosy90_rotenone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/cosy90_rotenone.m`
- Signature: `cosy90_rotenone()`
- Total lines: 84

## Purpose

COSY-90 spectrum of rotenone. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- COSY-90 spectrum of rotenone.
- Calculation time: minutes
- Spin system
- Interactions
- Algorithmic options
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
