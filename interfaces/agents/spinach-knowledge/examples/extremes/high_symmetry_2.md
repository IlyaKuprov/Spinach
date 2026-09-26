# examples/extremes/high_symmetry_2.m

- Signature: `high_symmetry_2()`

## Purpose

31P NMR spectrum of a large and highly symmetric spin system with two tert-butyl groups supplied by Eberhard Matern. Done by brute force time propagation in Hilbert space. WARNING: needs 32+ CPU cores and 128+ GB of RAM. Run time on the above: hours

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 31P NMR spectrum of a large and highly symmetric spin system
- with two tert-butyl groups supplied by Eberhard Matern. Done
- by brute force time propagation in Hilbert space.
- WARNING: needs 32+ CPU cores and 128+ GB of RAM.
- Run time on the above: hours
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
