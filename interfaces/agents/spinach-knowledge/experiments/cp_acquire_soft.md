# experiments/cp_acquire_soft.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/cp_acquire_soft.m`
- Signature: `fid=cp_acquire_soft(spin_system,parameters,H,R,K)`
- Total lines: 159

## Purpose

Cross-polarisation experiment in the rotating frame, followed by time-domain FID acquisition. The CP stage is preceded by wiping of the low-gamma spins and followed by FID acquisition with deco- upling of the high-gamma spins. Syntax: fid=cp_acquire_soft(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 54-55: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 57-58: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 60-62: Wipe the state of 13C; implemented by `[~,rho]=decouple(spin_system,[],parameters.rho0, parameters.spins(2))`.
- Lines 64-65: Build and project 1H and 13C control operators; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 72-74: Apply the 90-degree pulse on 1H along +X; implemented by `rho=step(spin_system,L+2*pi*parameters.hi_pwr*Hx, rho,1/(4*parameters.hi_pwr))`.
- Lines 76-80: Run the CP contact time evolution: irradiation of 1H along -Y, and of 13C along +X; implemented by `rho=evolution(spin_system,L-2*pi*parameters.cp_pwr(1)*Hy +2*pi*parameters.cp_pwr(2)*Cx, [],rho,parameters.cp_dur,1,'final')`.
- Lines 82-83: Wipe the state of 1H and apply 1H decoupling; implemented by `[L,rho]=decouple(spin_system,L,rho,parameters.spins(1))`.
- Lines 85-87: Run the acquisition; implemented by `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Key state/data transformations

- Lines 58: computes `L` using `L=H+1i*R+1i*K`.
- Lines 61-62: computes `[~,rho]` using `[~,rho]=decouple(spin_system,[],parameters.rho0, parameters.spins(2))`.
- Lines 65: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 66: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 67: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 73-74: computes `rho` using `rho=step(spin_system,L+2*pi*parameters.hi_pwr*Hx, rho,1/(4*parameters.hi_pwr))`.
- Lines 83: computes `[L,rho]` using `[L,rho]=decouple(spin_system,L,rho,parameters.spins(1))`.
- Lines 86-87: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 92: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.spins -working spins, a cell array of
- strings with high-gamma spin fi-
- rst, and low-gamma spin last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- channel, Hz
- parameters.cp_pwr -nutation frequencies on the two
- channels during the CP contact
- time, a two-element vector, Hz
- parameters.cp_dur -duration of the contact time, s
- parameters.rho0 -initial state, the state of the
- low-gamma spins will be wiped
- parameters.coil -detection state
- parameters.sweep -sweep width for the FID, Hz
- parameters.npoints -number of points in the FID
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- fid -signal detected on the coil state during
- system evolution

## Implementation structure

- Cross-polarisation experiment in the rotating frame, followed by
- time-domain FID acquisition. The CP stage is preceded by wiping
- of the low-gamma spins and followed by FID acquisition with deco-
- upling of the high-gamma spins. Syntax:
- fid=cp_acquire_soft(spin_system,parameters,H,R,K)
- parameters.spins -working spins, a cell array of
- strings with high-gamma spin fi-
- rst, and low-gamma spin last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- channel, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `iscell()`, `cellfun()`, `isscalar()`, `any()`, `ismember()`.
