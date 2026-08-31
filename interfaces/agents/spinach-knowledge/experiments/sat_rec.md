# experiments/sat_rec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/sat_rec.m`
- Signature: `fids=sat_rec(spin_system,parameters,H,R,K)`
- Total lines: 130

## Purpose

Saturation-recovery pulse sequence with analytical saturation (just the unit state as the initial condition). Syntax: fids=sat_rec(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 42-43: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 45-46: Get the initial condition; implemented by `rho0=unit_state(spin_system)`.
- Lines 48-49: Get the detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1})`.
- Lines 51-52: Get the pulse operator; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1}); Ly=(Lp-Lp')/2i`.
- Lines 54-57: Run relaxation periods; implemented by `rho_stack=evolution(spin_system,L,[],rho0, parameters.max_delay/parameters.n_delays, parameters.n_delays,'trajectory')`.
- Lines 59-60: Run a 90-degree pulse; implemented by `rho_stack=step(spin_system,Ly,rho_stack,pi/2)`.
- Lines 62-64: Run the detection period; implemented by `fids=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep, parameters.npoints-1,'observable')`.

### Key state/data transformations

- Lines 43: computes `L` using `L=H+1i*R+1i*K`.
- Lines 46: computes `rho0` using `rho0=unit_state(spin_system)`.
- Lines 49: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1})`.
- Lines 52: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1}); Ly=(Lp-Lp')/2i`.
- Lines 55-57: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho0, parameters.max_delay/parameters.n_delays, parameters.n_delays,'trajectory')`.
- Lines 63-64: computes `fids` using `fids=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep, parameters.npoints-1,'observable')`.

### Local helper functions

- Line 69: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.sweep spectrum sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.max_delay longest relaxation delay
- parameters.n_delays number of relaxation delays to run
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fids -free induction decays for each delay starting from zero,
- a matrix with individual FIDs in columns
- Note: the relaxation superoperator must be thermalised.
- Zak El-Machachi

## Implementation structure

- Saturation-recovery pulse sequence with analytical saturation (just
- the unit state as the initial condition). Syntax:
- fids=sat_rec(spin_system,parameters,H,R,K)
- parameters.sweep spectrum sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.max_delay longest relaxation delay
- parameters.n_delays number of relaxation delays to run
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `unit_state()`, `state()`, `operator()`, `evolution()`, `step()`, `ismatrix()`, `all()`, `isfield()`, `iscell()`, `ischar()`, `ismember()`, `isscalar()`, `strcmp()`.
