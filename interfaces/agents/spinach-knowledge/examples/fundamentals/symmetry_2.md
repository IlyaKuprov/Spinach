# examples/fundamentals/symmetry_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/symmetry_2.m`
- Signature: `symmetry_2()`
- Total lines: 68

## Purpose

Pulse-acquire NMR spectrum of a highly symmetric spin system provided by Andres Castillo. Uses the fully sym- metric irreducible representation of S3(x)S3(x)S3 per- mutation symmetry group.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Pulse-acquire NMR spectrum of a highly symmetric spin
- system provided by Andres Castillo. Uses the fully sym-
- metric irreducible representation of S3(x)S3(x)S3 per-
- mutation symmetry group.
- Spin system specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
