# examples/fundamentals/symmetry_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/symmetry_3.m`
- Signature: `symmetry_3()`
- Total lines: 57

## Purpose

1H NMR spectrum of valine. Uses the fully symmetric irreducible representation of the S3(x)S3 group.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- 1H NMR spectrum of valine. Uses the fully symmetric
- irreducible representation of the S3(x)S3 group.
- System specificaton
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
