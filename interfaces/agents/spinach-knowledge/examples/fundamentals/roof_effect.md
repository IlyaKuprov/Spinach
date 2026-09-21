# examples/fundamentals/roof_effect.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/roof_effect.m`
- Signature: `roof_effect()`
- Total lines: 67

## Purpose

Roof effect in a strongly J-coupled two-spin system.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Roof effect in a strongly J-coupled two-spin system.
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
- Get the figure going
- Loop over line positions
- Update the Zeeman frequencies
- Run the simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `kfigure()`, `scale_figure()`, `spin()`, `liquid()`, `apodisation()`, `fftshift()`, `subplot()`, `plot_1d()`, `kylabel()`.
