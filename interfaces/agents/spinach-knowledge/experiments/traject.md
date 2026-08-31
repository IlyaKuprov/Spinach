# experiments/traject.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/traject.m`
- Signature: `traj=traject(spin_system,parameters,H,R,K)`
- Total lines: 97

## Purpose

Simple forward time evolution trajectory. Syntax: traj=traject(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 35-36: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 38-39: Apply the decoupling; implemented by `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 41-43: Run the evolution and watch the coil state; implemented by `traj=evolution(spin_system,L,[],parameters.rho0,1/parameters.sweep, parameters.npoints-1,'trajectory')`.

### Key state/data transformations

- Lines 36: computes `L` using `L=H+1i*R+1i*K`.
- Lines 39: computes `[L,parameters.rho0]` using `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 42-43: computes `traj` using `traj=evolution(spin_system,L,[],parameters.rho0,1/parameters.sweep, parameters.npoints-1,'trajectory')`.

### Local helper functions

- Line 48: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the trajectory
- parameters.rho0 initial state
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- traj -system trajectory, a bookshelf stack of state vectors

## Implementation structure

- Simple forward time evolution trajectory. Syntax:
- traj=traject(spin_system,parameters,H,R,K)
- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the trajectory
- parameters.rho0 initial state
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- traj -system trajectory, a bookshelf stack of state vectors
- Check consistency
- Compose Liouvillian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `spins()`, `iscell()`, `any()`, `cellfun()`, `ismember()`.
