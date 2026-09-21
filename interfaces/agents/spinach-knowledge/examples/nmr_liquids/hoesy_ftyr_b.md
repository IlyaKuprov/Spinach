# examples/nmr_liquids/hoesy_ftyr_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hoesy_ftyr_b.m`
- Signature: `hoesy_ftyr_b()`
- Total lines: 74

## Purpose

(19F) -> (1H) HOESY spectrum of fluorotyrosine. This is not the right way to run this sequence in proteins because aro- matic 19F T2 is short, but 19F is being phase-encoded. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- (19F) -> (1H) HOESY spectrum of fluorotyrosine. This is not
- the right way to run this sequence in proteins because aro-
- matic 19F T2 is short, but 19F is being phase-encoded.
- Calculation time: minutes
- Read 3-fluorotyrosine DFT calculation
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
