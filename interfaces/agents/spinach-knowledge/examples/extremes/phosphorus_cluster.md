# examples/extremes/phosphorus_cluster.m

- Signature: `phosphorus_cluster()`

## Purpose

Phosphorus system simulation for Gerhard Hagele. Done by brute force Liouville space time propagation. WARNING: needs 32+ CPU cores, 128+ GB of RAM and a strong FP64 capable Nvidia GPU. Run time on the above: hours

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Phosphorus system simulation for Gerhard Hagele. Done by brute
- force Liouville space time propagation.
- WARNING: needs 32+ CPU cores, 128+ GB of RAM and
- a strong FP64 capable Nvidia GPU.
- Run time on the above: hours
- Magnet induction
- Isotopes
- Chemical shifts
- J-couplings
- Basis set
- Symmetry
- Greedy parallelisation
