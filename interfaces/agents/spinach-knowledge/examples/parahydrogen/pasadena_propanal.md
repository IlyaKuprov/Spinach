# examples/parahydrogen/pasadena_propanal.m

- Signature: `pasadena_propanal()`

## Purpose

PASADENA experiment simulation for the parahydrogenation of acrolein into propanal. Calculation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- PASADENA experiment simulation for the parahydrogenation of acrolein
- into propanal.
- Calculation time: seconds
- Spin system
- Magnetic field
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
