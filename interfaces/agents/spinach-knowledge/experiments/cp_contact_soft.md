# experiments/cp_contact_soft.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/cp_contact_soft.m`
- Signature: `contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)`
- Total lines: 150

## Purpose

Cross-polarisation experiment in the rotating frame. Applies a soft pi/2 pulse using the specified operators, then evolves the system with the specified spin-lock terms added to the Liovilli- an. The contact curve is returned. Syntax: contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 54-55: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 57-58: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 60-62: Wipe the state of 13C; implemented by `[~,rho]=decouple(spin_system,[],parameters.rho0, parameters.spins(2))`.
- Lines 64-65: Build and project 1H and 13C control operators; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 72-74: Apply the 90-degree pulse on 1H along +X; implemented by `rho=step(spin_system,L+2*pi*parameters.hi_pwr*Hx, rho,1/(4*parameters.hi_pwr))`.
- Lines 76-81: Run the CP contact time evolution: irradiation of 1H along -Y, and of 13C along +X; implemented by `contact_curve=evolution(spin_system,L-2*pi*parameters.cp_pwr(1)*Hy +2*pi*parameters.cp_pwr(2)*Cx, parameters.coil,rho,parameters.timestep, parameters.nsteps,'observable')`.

### Key state/data transformations

- Lines 58: computes `L` using `L=H+1i*R+1i*K`.
- Lines 61-62: computes `[~,rho]` using `[~,rho]=decouple(spin_system,[],parameters.rho0, parameters.spins(2))`.
- Lines 65: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 66: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 67: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 73-74: computes `rho` using `rho=step(spin_system,L+2*pi*parameters.hi_pwr*Hx, rho,1/(4*parameters.hi_pwr))`.
- Lines 78-81: computes `contact_curve` using `contact_curve=evolution(spin_system,L-2*pi*parameters.cp_pwr(1)*Hy +2*pi*parameters.cp_pwr(2)*Cx, parameters.coil,rho,parameters.timestep, parameters.nsteps,'observable')`.

### Local helper functions

- Line 86: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.spins -working spins, a cell array of
- strings with high-gamma spins
- first and low-gamma spins last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- spins, Hz
- parameters.cp_pwr -nutation frequencies on the two
- channels during the CP contact
- time, a two-element vector, Hz
- parameters.timestep -time step of the CP contact ti-
- me, seconds
- parameters.nsteps -number of time steps to take
- during the CP contact time
- parameters.rho0 -initial state, the state of the
- low-gamma spins will be wiped
- parameters.coil -detection state vector
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- contact_curve -contact curve detected on the coil
- state specified in parameters.coil

## Implementation structure

- Cross-polarisation experiment in the rotating frame. Applies a
- soft pi/2 pulse using the specified operators, then evolves the
- system with the specified spin-lock terms added to the Liovilli-
- an. The contact curve is returned. Syntax:
- contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)
- parameters.spins -working spins, a cell array of
- strings with high-gamma spins
- first and low-gamma spins last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- spins, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `iscell()`, `cellfun()`, `isscalar()`, `isrow()`, `any()`.
