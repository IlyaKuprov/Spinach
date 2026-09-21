# examples/nmr_liquids/pansy_triple_ch.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pansy_triple_ch.m`
- Signature: `pansy_triple_ch()`
- Total lines: 111

## Purpose

Triple-channel PANSY experiment on glycine with natural content of 13C isotope. Coordinates, shieldings, and J- couplings computed with DFT. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Triple-channel PANSY experiment on glycine with natural
- content of 13C isotope. Coordinates, shieldings, and J-
- couplings computed with DFT.
- Calculation time: seconds
- Read the spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Spinach housekeeping
- Sequence parameters
- Preallocate the answer
- Generate isotopomers
- Loop over isotopomers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_offset()`, `plot_sweep()`, `plot_zerofill()`, `plot_spins()`.
