# examples/kinetics/exchange_symmetric.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/exchange_symmetric.m`
- Signature: `exchange_symmetric()`
- Total lines: 51

## Purpose

Two-spin symmetric chemical exchange pattern. Calculation time: seconds.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Two-spin symmetric chemical exchange pattern.
- Calculation time: seconds.
- System specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
