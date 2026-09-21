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
