# examples/giant_spin/quartet_magn.m

- Signature: `quartet_magn()`

## Purpose

Sample magnetisation during a finite-speed magnetic field sweep for a spin-3/2 particle with a zero-field splitting. Calculation time: seconds

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Sample magnetisation during a finite-speed magnetic field
- sweep for a spin-3/2 particle with a zero-field splitting.
- Calculation time: seconds
- This must be set to 1 Tesla
- Particle
- Zeeman tensor
- Zero-field splitting
- Formalism and basis set
- Temperature
- Spinach housekeeping
- Experiment parameters
- Run the field scan
