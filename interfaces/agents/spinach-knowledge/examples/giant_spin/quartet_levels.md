# examples/giant_spin/quartet_levels.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/quartet_levels.m`
- Signature: `quartet_levels()`
- Total lines: 42

## Purpose

Energy levels magnetic field scan for a spin-3/2 particle with a zero-field splitting. Calculation time: seconds

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Energy levels magnetic field scan for a spin-3/2 particle
- with a zero-field splitting.
- Calculation time: seconds
- This must be set to 1 Tesla
- Particle
- Zeeman tensor
- Zero-field splitting
- Formalism and basis set
- Spinach housekeeping
- Experiment parameters
- Run the field scan

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `zfs2mat()`, `create()`, `basis()`, `fieldscan_enlev()`.
