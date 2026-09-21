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
