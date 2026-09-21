# examples/nmr_liquids/mqs_propanol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/mqs_propanol.m`
- Signature: `mqs_propanol()`
- Total lines: 101

## Purpose

Multiple-quantum NMR experiment for a propanol spin system. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Multiple-quantum NMR experiment for a propanol
- spin system.
- Calculation time: seconds
- Magnet field
- Chemical shifts
- 2J couplings
- 3J couplings
- Basis set
- Algorithmic options
- Spinach housekeeping
- Initial and detection states
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `kfigure()`, `scale_figure()`, `tau_max()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `subplot()`, `plot_2d()`, `num2str()`.
