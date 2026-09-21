# examples/nmr_liquids/hmbc_camphor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hmbc_camphor.m`
- Signature: `hmbc_camphor()`
- Total lines: 74

## Purpose

HMBC spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings computed witt DFT. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- HMBC spectrum of camphor with natural content of 13C isotope.
- Coordinates, shielding anisotropies and J-couplings computed
- witt DFT.
- Calculation time: seconds
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Algorithmic options
- Basis set
- Sequence parameters
- Create the spin system structure
- Generate isotopomers
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
