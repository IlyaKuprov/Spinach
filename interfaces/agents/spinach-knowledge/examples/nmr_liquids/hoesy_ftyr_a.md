# examples/nmr_liquids/hoesy_ftyr_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hoesy_ftyr_a.m`
- Signature: `hoesy_ftyr_a()`
- Total lines: 75

## Purpose

(1H) -> (19F) HOESY spectrum of fluorotyrosine, with the magneti- sation transfer direction picked so as to minimise the time that 19F spends in the transverse plane. This is the only way to run this sequence in proteins because aromatic 19F T2 is short. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- (1H) -> (19F) HOESY spectrum of fluorotyrosine, with the magneti-
- sation transfer direction picked so as to minimise the time that
- 19F spends in the transverse plane. This is the only way to run
- this sequence in proteins because aromatic 19F T2 is short.
- Calculation time: minutes
- Read 3-fluorotyrosine DFT calculation
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
