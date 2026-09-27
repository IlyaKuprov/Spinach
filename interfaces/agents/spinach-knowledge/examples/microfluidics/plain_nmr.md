# examples/microfluidics/plain_nmr.m

- Signature: `plain_nmr()`

## Purpose

NMR spectrum of the reaction mixture in the absence of chemical kinetics and spatial dynamics.

## Physical / mathematical content

- This is a homogeneous liquid-state NMR calculation, not a spatial microfluidics simulation. It imports the Diels–Alder reaction spin system, sets equal concentrations for the four chemical species with no solvent, and acquires a proton spectrum without chemical kinetics or spatial dynamics.
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
