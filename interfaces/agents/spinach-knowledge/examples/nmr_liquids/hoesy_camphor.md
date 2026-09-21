# examples/nmr_liquids/hoesy_camphor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hoesy_camphor.m`
- Signature: `hoesy_camphor()`
- Total lines: 89

## Purpose

13C{1H} HOESY spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings computed with DFT. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- 13C{1H} HOESY spectrum of camphor with natural content of 13C isotope.
- Coordinates, shielding anisotropies and J-couplings computed with DFT.
- Calculation time: minutes
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Generate isotopomers
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
