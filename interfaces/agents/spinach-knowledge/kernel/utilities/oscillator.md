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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(parameters)`.
- Lines 37-38: Second derivative operator; implemented by `d2_dx2=((parameters.n_points-1)/parameters.box_size)^2*fdmat(parameters.n_points,5,2)`.
- Lines 40-41: Coordinate grid; implemented by `xgrid=linspace(-parameters.box_size/2,parameters.box_size/2,parameters.n_points)'`.
- Lines 43-44: Coordinate operator; implemented by `X_oscl=spdiags(xgrid,0,parameters.n_points,parameters.n_points)`.
- Lines 46-47: Boundary conditions; implemented by `d2_dx2(:,1)=0; d2_dx2(:,end)=0; d2_dx2(1,:)=0; d2_dx2(end,:)=0`.
- Lines 49-52: Hamiltonian; implemented by `H_oscl=-(1/(2*parameters.par_mass))*d2_dx2+ (parameters.frc_cnst/2)*X_oscl^2+ parameters.par_mass*parameters.grv_cnst*X_oscl`.

### Key state/data transformations

- Lines 38: computes `d2_dx2` using `d2_dx2=((parameters.n_points-1)/parameters.box_size)^2*fdmat(parameters.n_points,5,2)`.
- Lines 41: computes `xgrid` using `xgrid=linspace(-parameters.box_size/2,parameters.box_size/2,parameters.n_points)'`.
- Lines 44: computes `X_oscl` using `X_oscl=spdiags(xgrid,0,parameters.n_points,parameters.n_points)`.
- Lines 47: computes `d2_dx2(:,1)` using `d2_dx2(:,1)=0; d2_dx2(:,end)=0; d2_dx2(1,:)=0; d2_dx2(end,:)=0`.
- Lines 50-52: computes `H_oscl` using `H_oscl=-(1/(2*parameters.par_mass))*d2_dx2+ (parameters.frc_cnst/2)*X_oscl^2+ parameters.par_mass*parameters.grv_cnst*X_oscl`.

### Local helper functions

- Line 57: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isfield(parameters,'frc_cnst')`.
  - Representative operation: `error('force constant must be specified in parameters.frc_cnst variable.')`.

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
