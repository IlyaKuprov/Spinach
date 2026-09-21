# examples/nmr_liquids/pa_styrene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_styrene.m`
- Signature: `pa_styrene()`
- Total lines: 64

## Purpose

1H NMR spectrum of styrene, the calculation is performed in Hilbert space to demonstrate parallel propagation described in: Calculation time: seconds.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H NMR spectrum of styrene, the calculation is performed in Hilbert
- space to demonstrate parallel propagation described in:
- Calculation time: seconds.
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
