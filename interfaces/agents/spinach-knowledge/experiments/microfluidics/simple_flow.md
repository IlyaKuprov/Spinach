# experiments/microfluidics/simple_flow.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/microfluidics/simple_flow.m`
- Signature: `traj=simple_flow(spin_system,parameters,H,R,K,~,F)`
- Total lines: 77

## Purpose

Simple forward evolution experiment for the microfluidics module; trajectory is returned. Syntax: traj=simple_flow(spin_system,parameters,H,R,K,~,F) This sequence must be called from the meshflow() context, which would provide H, R, K, G, and F. Because gradients are not being used, the G input is ignored.

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(parameters,H,R,K,F)`.
- Lines 37-38: Compose Liouvillian; implemented by `L=H+1i*F+1i*R+1i*K`.
- Lines 40-42: Run the evolution and watch the coil state; implemented by `traj=evolution(spin_system,L,[],parameters.rho0,parameters.dt, parameters.npoints-1,'trajectory')`.

### Key state/data transformations

- Lines 38: computes `L` using `L=H+1i*F+1i*R+1i*K`.
- Lines 41-42: computes `traj` using `traj=evolution(spin_system,L,[],parameters.rho0,parameters.dt, parameters.npoints-1,'trajectory')`.

### Local helper functions

- Line 47: `grumble()` — `function grumble(parameters,H,R,K,F)`.
  - Representative operation: `if ~isfield(parameters,'rho0')`.
  - Representative operation: `error('initial state must be specified in parameters.rho0 variable.')`.

## Parameters / inputs

- parameters.npoints -number of points in
- the trajectory
- parameters.rho0 -initial state in Fokker-
- Planck space
- parameters.dt -trajectory time step

## Outputs

- traj -trajectory in the Fokker-Planck space
- Notes: to convert Fokker-Planck space trajectory into R3
- or Liouville space, use fpl2phan and fpl2rho func-
- tions.

## Implementation structure

- Simple forward evolution experiment for the microfluidics
- module; trajectory is returned. Syntax:
- traj=simple_flow(spin_system,parameters,H,R,K,~,F)
- This sequence must be called from the meshflow() context,
- which would provide H, R, K, G, and F. Because gradients
- are not being used, the G input is ignored.
- parameters.npoints -number of points in
- the trajectory
- parameters.rho0 -initial state in Fokker-
- Planck space
- parameters.dt -trajectory time step
- traj -trajectory in the Fokker-Planck space

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `evolution()`, `isfield()`, `ismatrix()`, `isscalar()`.
