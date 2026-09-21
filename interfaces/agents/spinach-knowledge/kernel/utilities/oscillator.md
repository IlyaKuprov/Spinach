# kernel/utilities/oscillator.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/oscillator.m`
- Signature: `[H_oscl,X_oscl,xgrid]=oscillator(parameters)`
- Total lines: 100

## Purpose

Harmonic oscillator infrastructure in 1D. Syntax: [H_oscl,X_oscl,xgrid]=oscillator(parameters)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.frc_cnst -force constant, N/m
- parameters.par_mass -particle mass, kg
- parameters.grv_cnst -gravitational acceleration. m/s^2
- parameters.n_points -number of discretization points
- parameters.box_size -oscillator box size, m

## Outputs

- H_oscl -oscillator Hamiltonian, Joules
- X_oscl -oscillator X operator, m
- xgrid -X coordinate grid, m
- Note: gravitation is directed along the X axis. Finite difference
- derivative operators are used.

## Implementation structure

- Harmonic oscillator infrastructure in 1D. Syntax:
- [H_oscl,X_oscl,xgrid]=oscillator(parameters)
- parameters.frc_cnst -force constant, N/m
- parameters.par_mass -particle mass, kg
- parameters.grv_cnst -gravitational acceleration. m/s^2
- parameters.n_points -number of discretization points
- parameters.box_size -oscillator box size, m
- H_oscl - oscillator Hamiltonian, Joules
- X_oscl - oscillator X operator, m
- xgrid - X coordinate grid, m
- Note: gravitation is directed along the X axis. Finite difference
- derivative operators are used.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fdmat()`, `spdiags()`, `d2_dx2()`, `isfield()`.
