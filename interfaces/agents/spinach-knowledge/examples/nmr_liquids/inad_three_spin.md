# examples/nmr_liquids/inad_three_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inad_three_spin.m`
- Signature: `inad_three_spin()`
- Total lines: 50

## Purpose

INADEQUATE spectrum of a three-spin system with J-coupling between two spins only. The sequence selects double-quantum coherence from coupled 13C pairs and converts it back for detection. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- INADEQUATE spectrum of a three-spin system with J-coupling between
- two spins only. The sequence selects double-quantum coherence
- from coupled 13C pairs and converts it back for detection.
- Calculation time: seconds
- Magnet field
- Spin system and interactions
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Processing
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
