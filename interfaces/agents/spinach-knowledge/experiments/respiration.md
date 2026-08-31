# experiments/respiration.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/respiration.m`
- Signature: `fid=respiration(spin_system,parameters,H,R,K)`
- Total lines: 154

## Purpose

RESPIRATION cross-polarisation method described in the paper from the Aarhus group (http://dx.doi.org/10.1021/jz3000905). Syntax: fid=respiration(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 47-48: Generate and project pulse operators; implemented by `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 54-55: Build the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 57-58: Initial condition; implemented by `rho=parameters.rho0`.
- Lines 60-61: RESPIRATION loop; implemented by `for n=1:parameters.nloops`.
- Lines 63-65: Update the user; implemented by `report(spin_system,['RESPIRATION loop ' num2str(n) '/' num2str(parameters.nloops) ' '])`.
- Lines 67-68: +X pulse on protons; implemented by `rho=step(spin_system,L+2*pi*2*parameters.rate*Hx,rho,1/(2*parameters.rate))`.
- Lines 70-71: -X pulse on protons; implemented by `rho=step(spin_system,L-2*pi*2*parameters.rate*Hx,rho,1/(2*parameters.rate))`.
- Lines 73-74: Ideal theta pulse on both; implemented by `rho=step(spin_system,Cx+Hx,rho,parameters.theta)`.
- Lines 78-79: Acquisition; implemented by `[L,rho]=decouple(spin_system,L,rho,parameters.spins(1))`.

### Control flow inferred from the code

- Line 61: `for` loop over `n=1:parameters.nloops`.

### Key state/data transformations

- Lines 48: computes `Hp` using `Hp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 49: computes `Cp` using `Cp=operator(spin_system,'L+',parameters.spins{2})`.
- Lines 52: computes `Hx` using `Hx=(Hp+Hp')/2; Cx=(Cp+Cp')/2`.
- Lines 55: computes `L` using `L=H+1i*R+1i*K`.
- Lines 58: computes `rho` using `rho=parameters.rho0`.
- Lines 79: computes `[L,rho]` using `[L,rho]=decouple(spin_system,L,rho,parameters.spins(1))`.
- Lines 80-81: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 86: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.nloops number of RESPIRATION loops
- parameters.theta the angle of the ideal pulse
- at the end of each loop
- parameters.rate RESPIRATION pulse train rate, Hz
- parameters.spins working spins, e.g. {'1H','13C'}
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay as seen by the state specified
- in parameters parameters.coil

## Implementation structure

- RESPIRATION cross-polarisation method described in the paper from
- the Aarhus group (http://dx.doi.org/10.1021/jz3000905). Syntax:
- fid=respiration(spin_system,parameters,H,R,K)
- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.nloops number of RESPIRATION loops
- parameters.theta the angle of the ideal pulse
- at the end of each loop
- parameters.rate RESPIRATION pulse train rate, Hz
- parameters.spins working spins, e.g. {'1H','13C'}

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `report()`, `num2str()`, `step()`, `decouple()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `isscalar()`, `iscell()`, `cellfun()`, `any()`, `ismember()`.
