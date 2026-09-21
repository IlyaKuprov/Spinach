# examples/nmr_liquids/pa_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_sucrose.m`
- Signature: `pa_sucrose()`
- Total lines: 61

## Purpose

1H NMR spectrum of sucrose (magnetic parameters read in from a DFT calculation), including Redfield relaxation superoperator. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H NMR spectrum of sucrose (magnetic parameters read in from a DFT
- calculation), including Redfield relaxation superoperator.
- Calculation time: seconds
- Read the spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
