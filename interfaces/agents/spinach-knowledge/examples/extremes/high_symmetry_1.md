# examples/extremes/high_symmetry_1.m

- Signature: `high_symmetry_1()`

## Purpose

1H NMR spectrum of a large and highly symmetric spin system with two tert-butyl groups supplied by Eberhard Matern. Done by brute force time-domain propagation in Hilbert space. WARNING: needs 32+ CPU cores and 128+ GB of RAM. Run time on the above: hours

## Physical / mathematical content

- The 1H NMR spectrum is calculated for a large, highly symmetric system containing two tert-butyl groups and phosphorus nuclei.
- The script calculates a time-domain NMR free-induction decay and Fourier-transforms it to obtain the plotted spectrum.

## Numerical / algorithmic content

- The script performs brute-force time-domain propagation in Hilbert space to calculate the FID and spectrum.

## Implementation structure

- 1H NMR spectrum of a large and highly symmetric spin system
- with two tert-butyl groups supplied by Eberhard Matern. Done
- by brute force time-domain propagation in Hilbert space.
- WARNING: needs 32+ CPU cores and 128+ GB of RAM.
- Run time on the above: hours
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
