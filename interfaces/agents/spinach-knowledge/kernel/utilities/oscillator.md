# kernel/utilities/oscillator.m

- Signature: `[H_oscl,X_oscl,xgrid]=oscillator(parameters)`

## Purpose

Harmonic oscillator infrastructure in 1D. Syntax: [H_oscl,X_oscl,xgrid]=oscillator(parameters)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

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
