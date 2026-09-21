# examples/nmr_liquids/cosy90_derome.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/cosy90_derome.m`
- Signature: `cosy90_derome()`
- Total lines: 59

## Purpose

Figure 8.26 from Andrew Derome's "Modern NMR Techniques for Chemistry Research". Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Figure 8.26 from Andrew Derome's "Modern NMR Techniques
- for Chemistry Research".
- Calculation time: seconds
- Magnet field
- Spin system and interactions
- Basis set
- Algorithmic options
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
