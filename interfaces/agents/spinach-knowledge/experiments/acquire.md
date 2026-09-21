# experiments/acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/acquire.m`
- Signature: `fid=acquire(spin_system,parameters,H,R,K)`
- Total lines: 158

## Purpose

Simple forward time evolution with signal acquisition. Syntax: fid=acquire(spin_system,parameters,H,R,K)

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
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- parameters.homodec_oper operator to add to the Liouvillian at
- the detection stage
- parameters.homodec_pwr power coefficient for the operator, Hz
- parameters.dead_time the system will be evolved for this
- time (seconds) before the signal
- acquisition begins
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay as seen by the state specified
- in parameters parameters.coil

## Implementation structure

- Simple forward time evolution with signal acquisition. Syntax:
- fid=acquire(spin_system,parameters,H,R,K)
- parameters.sweep sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.decouple spins to decouple, e.g. {'15N','13C'}
- parameters.homodec_oper operator to add to the Liouvillian at
- the detection stage
- parameters.homodec_pwr power coefficient for the operator, Hz
- parameters.dead_time the system will be evolved for this
- time (seconds) before the signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `isfield()`, `speye()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isscalar()`, `spins()`, `iscell()`, `any()`, `cellfun()`, `ismember()`.
