# examples/nmr_solids/mqmas_nqi.m

- Signature: `mqmas_nqi()`

## Purpose

Rotor-synchronous MQMAS spectrum of a 87Rb compound, transmitter set to the isotropic chemical shift. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Rotor-synchronous MQMAS spectrum of a 87Rb compound,
- transmitter set to the isotropic chemical shift.
- Calculation time: minutes
- System specification: just the NQI
- Formalism and basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform
- Plotting
