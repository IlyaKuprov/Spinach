# examples/nmr_liquids/clip_hsqc_camphor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/clip_hsqc_camphor.m`
- Signature: `clip_hsqc_camphor()`
- Total lines: 87

## Purpose

CLIP-HSQC spectrum of camphor with natural content of 13C isotope. Coordinates, shielding anisotropies and J-couplings computed with DFT, isotropic chemical shifts taken from experimental data. Calculation time: minutes.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Implementation structure

- CLIP-HSQC spectrum of camphor with natural content of 13C isotope.
- Coordinates, shielding anisotropies and J-couplings computed with
- DFT, isotropic chemical shifts taken from experimental data.
- Calculation time: minutes.
- Spin system properties (vacuum DFT calculation)
- Magnet field
- Isotropic shift components come from the experiment
- Basis set
- Algorithmic options
- Sequence parameters
- Create the spin system structure
- Generate isotopomers

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `conj()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
