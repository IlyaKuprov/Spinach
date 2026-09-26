# examples/microfluidics/plain_nmr.m

- Signature: `plain_nmr()`

## Purpose

NMR spectrum of the reaction mixture in the absence of chemical kinetics and spatial dynamics.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- NMR spectrum of the reaction mixture in the absence of
- chemical kinetics and spatial dynamics.
- Import Diels-Alder cycloaddition
- Equal concentrations, no solvent
- Magnet field
- Greedy parallelisation
- Spinach housekeeping
- Sequence parameters -1H
- Simulation
- Apodisation
- Fourier transform
- Plotting
