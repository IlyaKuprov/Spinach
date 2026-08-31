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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 50-51: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 53-54: Compose the Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 56-57: Project pulse operator; implemented by `parameters.pulse_op=kron(speye(parameters.spc_dim),parameters.pulse_op)`.
- Lines 59-60: Apply the pulse; implemented by `rho=step(spin_system,parameters.pulse_op,parameters.rho0,parameters.pulse_angle)`.
- Lines 62-63: Run the echo stage if specified; implemented by `if isfield(parameters,'echo_time')`.
- Lines 65-66: Run the first evolution time in the echo train; implemented by `rho=evolution(spin_system,L,[],rho,parameters.echo_time,1,'final')`.
- Lines 68-69: Project pulse operator; implemented by `parameters.echo_oper=kron(speye(parameters.spc_dim),parameters.echo_oper)`.
- Lines 71-72: Run the echo pulse; implemented by `rho=step(spin_system,parameters.echo_oper,rho,parameters.echo_angle)`.
- Lines 74-75: Run the second evolution time in the echo train; implemented by `rho=evolution(spin_system,L,[],rho,parameters.echo_time,1,'final')`.
- Lines 79-80: Apply the decoupling; implemented by `[L,rho]=decouple(spin_system,L,rho,parameters.decouple)`.
- Lines 82-84: Run the evolution and watch the coil state; implemented by `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Control flow inferred from the code

- Line 63: conditional branch on `isfield(parameters,'echo_time')`.

### Key state/data transformations

- Lines 54: computes `L` using `L=H+1i*R+1i*K`.
- Lines 57: computes `parameters.pulse_op` using `parameters.pulse_op=kron(speye(parameters.spc_dim),parameters.pulse_op)`.
- Lines 60: computes `rho` using `rho=step(spin_system,parameters.pulse_op,parameters.rho0,parameters.pulse_angle)`.
- Lines 69: computes `parameters.echo_oper` using `parameters.echo_oper=kron(speye(parameters.spc_dim),parameters.echo_oper)`.
- Lines 80: computes `[L,rho]` using `[L,rho]=decouple(spin_system,L,rho,parameters.decouple)`.
- Lines 83-84: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 89: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `error('H, R and K arguments must be matrices.')`.

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
