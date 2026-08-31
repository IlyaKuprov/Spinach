# experiments/cpmg.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/cpmg.m`
- Signature: `fid=cpmg(spin_system,parameters,H,R,K)`
- Total lines: 118

## Purpose

CPMG echo train with detection. Syntax: fid=cpmg(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 32-33: Project the operator; implemented by `parameters.pulse_op=kron(speye(parameters.spc_dim),parameters.pulse_op)`.
- Lines 35-36: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 38-39: Excitation pulse; implemented by `rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/2)`.
- Lines 41-43: Run a half-echo; implemented by `traj=evolution(spin_system,L,[],rho, parameters.timestep,parameters.npoints-1,'trajectory')`.
- Lines 45-46: Record the half-echo; implemented by `fid=parameters.coil'*traj`.
- Lines 48-49: CPMG loop; implemented by `for n=1:parameters.nloops`.
- Lines 51-52: Apply the pulse; implemented by `traj(:,end)=step(spin_system,parameters.pulse_op,traj(:,end),pi)`.
- Lines 54-56: Run the echo; implemented by `traj=evolution(spin_system,L,[],traj(:,end), parameters.timestep,2*parameters.npoints-1,'trajectory')`.
- Lines 58-59: Record the echo; implemented by `fid=[fid parameters.coil'*traj]; %#ok<AGROW>`.

### Control flow inferred from the code

- Line 49: `for` loop over `n=1:parameters.nloops`.

### Key state/data transformations

- Lines 33: computes `parameters.pulse_op` using `parameters.pulse_op=kron(speye(parameters.spc_dim),parameters.pulse_op)`.
- Lines 36: computes `L` using `L=H+1i*R+1i*K`.
- Lines 39: computes `rho` using `rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/2)`.
- Lines 42-43: computes `traj` using `traj=evolution(spin_system,L,[],rho, parameters.timestep,parameters.npoints-1,'trajectory')`.
- Lines 46: computes `fid` using `fid=parameters.coil'*traj`.
- Lines 52: computes `traj(:,end)` using `traj(:,end)=step(spin_system,parameters.pulse_op,traj(:,end),pi)`.

### Local helper functions

- Line 66: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `error('H, R and K arguments must be matrices.')`.

## Parameters / inputs

- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_op -pulse operator
- parameters.nloops -number of CPMG loops
- parameters.timestep -time step
- parameters.npoints -number of steps per half-echo

## Outputs

- fid -free induction decay throughout the sequence

## Implementation structure

- CPMG echo train with detection. Syntax:
- fid=cpmg(spin_system,parameters,H,R,K)
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_op -pulse operator
- parameters.nloops -number of CPMG loops
- parameters.timestep -time step
- parameters.npoints -number of steps per half-echo
- fid -free induction decay throughout the sequence
- Check consistency
- Project the operator
- Compose Liouvillian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `speye()`, `step()`, `evolution()`, `traj()`, `ismatrix()`, `all()`, `isfield()`.
