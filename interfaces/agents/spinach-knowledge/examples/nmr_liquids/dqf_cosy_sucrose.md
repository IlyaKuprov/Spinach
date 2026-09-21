# examples/nmr_liquids/dqf_cosy_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/dqf_cosy_sucrose.m`
- Signature: `dqf_cosy_sucrose()`
- Total lines: 62

## Purpose

DQF-COSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- DQF-COSY spectrum of sucrose (magnetic parameters computed with DFT).
- Calculation time: minutes
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodization
- F2 Fourier transform
- Form States signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
