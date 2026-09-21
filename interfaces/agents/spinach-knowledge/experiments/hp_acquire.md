# experiments/hp_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/hp_acquire.m`
- Signature: `fid=hp_acquire(spin_system,parameters,H,R,K)`
- Total lines: 167

## Purpose

Standard pulse-acquire sequence with a hard pulse. The user must sup- ply the pulse operator, the pulse duration and the initial condition. Echo detection may optionally be used. Syntax: fid=hp_acquire(spin_system,parameters,H,R,K)

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
- parameters.pulse_op pulse operator
- parameters.pulse_angle pulse angle, radians
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- (sphten-liouv formalism only)
- parameters.echo_time optional echo time for echo detection
- (echo_time -pulse -echo_time -fid)
- parameters.echo_oper optional pulse operator for echo
- detection
- parameters.echo_angle optional pulse angle for echo detection
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay as seen by the state specified
- in parameters parameters.coil

## Implementation structure

- Standard pulse-acquire sequence with a hard pulse. The user must sup-
- ply the pulse operator, the pulse duration and the initial condition.
- Echo detection may optionally be used. Syntax:
- fid=hp_acquire(spin_system,parameters,H,R,K)
- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.pulse_op pulse operator
- parameters.pulse_angle pulse angle, radians
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- (sphten-liouv formalism only)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `step()`, `isfield()`, `evolution()`, `decouple()`, `ismatrix()`, `all()`, `spins()`, `iscell()`, `any()`, `cellfun()`, `ismember()`, `isscalar()`.
