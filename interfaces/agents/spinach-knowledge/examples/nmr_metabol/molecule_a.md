# examples/nmr_metabol/molecule_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_metabol/molecule_a.m`
- Signature: `molecule_a()`
- Total lines: 47

## Purpose

1H NMR spectrum of a molecule from the GISSMO database. Calculation time: seconds

## Physical / mathematical content

- Metabolomics NMR examples. These files apply liquid-state NMR simulation workflows to small-molecule mixtures, spectral assignment, concentration inference, and database-style metabolite spin-system definitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H NMR spectrum of a molecule from the GISSMO database.
- Calculation time: seconds
- Import GISSMO dataset
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gissmo2spinach()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
