# examples/parahydrogen/pasadena_ethylbenzene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/pasadena_ethylbenzene.m`
- Signature: `pasadena_ethylbenzene()`
- Total lines: 70

## Purpose

PASADENA experiment simulation for the parahydrogenation of styrene into ethylbenzene. Set to reproduce the top trace of Fig 5 in Calculation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- PASADENA experiment simulation for the parahydrogenation of styrene
- into ethylbenzene. Set to reproduce the top trace of Fig 5 in
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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
