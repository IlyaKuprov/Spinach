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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 47-48: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 50-51: Apply the decoupling; implemented by `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 53-54: Apply the homodecoupling; implemented by `if isfield(parameters,'homodec_oper')`.
- Lines 56-57: Project into FP space; implemented by `parameters.homodec_oper=kron(speye(parameters.spc_dim),parameters.homodec_oper)`.
- Lines 59-60: Add to the Liouvillian; implemented by `L=L+2*pi*parameters.homodec_pwr*parameters.homodec_oper`.
- Lines 64-65: Run the dead time; implemented by `if isfield(parameters,'dead_time')`.
- Lines 69-71: Run the evolution and watch the coil state; implemented by `fid=evolution(spin_system,L,parameters.coil,parameters.rho0, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Control flow inferred from the code

- Line 54: conditional branch on `isfield(parameters,'homodec_oper')`.
- Line 65: conditional branch on `isfield(parameters,'dead_time')`.

### Key state/data transformations

- Lines 48: computes `L` using `L=H+1i*R+1i*K`.
- Lines 51: computes `[L,parameters.rho0]` using `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 57: computes `parameters.homodec_oper` using `parameters.homodec_oper=kron(speye(parameters.spc_dim),parameters.homodec_oper)`.
- Lines 66: computes `parameters.rho0` using `parameters.rho0=step(spin_system,L,parameters.rho0,parameters.dead_time)`.
- Lines 70-71: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,parameters.rho0, 1/parameters.sweep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 76: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

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
