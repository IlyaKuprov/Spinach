# examples/nmr_liquids/pa_rotenone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_rotenone.m`
- Signature: `pa_rotenone()`
- Total lines: 91

## Purpose

1H NMR spectrum of rotenone using T1/T2 relaxation model, magnetic parameters from: Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H NMR spectrum of rotenone using T1/T2 relaxation model,
- magnetic parameters from:
- Calculation time: seconds
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Relaxation model
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
