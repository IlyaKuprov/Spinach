# examples/extremes/fluoroisooctane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/fluoroisooctane.m`
- Signature: `fluoroisooctane()`
- Total lines: 96

## Purpose

A deliberately adversarial example from Art Bochevarov at Schodinger Inc. In this case, IK-2 approximation in Liou- ville space generates an exceedingly large basis set; the calculation must instead be performed in Hilbert space with permutation symmetry factorisation. Calculation time: hours.

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- A deliberately adversarial example from Art Bochevarov at
- Schodinger Inc. In this case, IK-2 approximation in Liou-
- ville space generates an exceedingly large basis set; the
- calculation must instead be performed in Hilbert space
- with permutation symmetry factorisation.
- Calculation time: hours.
- Magnet induction
- Isotopes
- Chemical shifts
- Larger J-couplings
- Smaller J-couplings, tert-butyl
- Smaller J-couplings, isopropyl

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `plot_1d()`.
