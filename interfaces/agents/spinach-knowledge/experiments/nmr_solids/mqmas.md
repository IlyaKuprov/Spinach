# experiments/nmr_solids/mqmas.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_solids/mqmas.m`
- Signature: `fid=mqmas(spin_system,parameters,H,R,K)`
- Total lines: 183

## Purpose

Rotor-synchronous MQMAS pulse sequence. Syntax: fid=mqmas(spin_system,parameters,H,R,K) This function should normally be called using singlerot.m context that would provide H, R, and K.

## Physical / mathematical content

- Solid-state pulse sequence implementations. The core ingredients are anisotropic Hamiltonians, rotor synchronisation, cross-polarisation, recoupling/decoupling, and powder or rotor-stack propagation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 46-47: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 49-50: Get pulse operators; implemented by `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 53-54: Apply the decoupling; implemented by `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 56-58: Run the first pulse; implemented by `rho=step(spin_system,L+parameters.pulse_amp(1)*Lx, parameters.rho0,parameters.pulse_dur(1))`.
- Lines 60-61: Do the coherence selection; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},parameters.mq_order}})`.
- Lines 63-65: Run the indirect dimension evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,1/parameters.rate, parameters.npoints(1)-1,'trajectory')`.
- Lines 67-69: Run the second pulse; implemented by `rho_stack=step(spin_system,L+parameters.pulse_amp(2)*Lx, rho_stack,parameters.pulse_dur(2))`.
- Lines 71-72: Do the coherence selection; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}})`.
- Lines 74-76: Run the direct dimension evolution; implemented by `fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.rate, parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 47: computes `L` using `L=H+1i*R+1i*K`.
- Lines 50: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 51: computes `Lx` using `Lx=(Lp+Lp')/2; Lx=kron(speye(parameters.spc_dim),Lx)`.
- Lines 54: computes `[L,parameters.rho0]` using `[L,parameters.rho0]=decouple(spin_system,L,parameters.rho0,parameters.decouple)`.
- Lines 57-58: computes `rho` using `rho=step(spin_system,L+parameters.pulse_amp(1)*Lx, parameters.rho0,parameters.pulse_dur(1))`.
- Lines 64-65: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,1/parameters.rate, parameters.npoints(1)-1,'trajectory')`.
- Lines 75-76: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.rate, parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))|| (~isequal(size(H),size(R),size(K)))|| (size(H,1)~=size(H,2))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))|| (~isequal(size(H),size(R),size(K)))|| (size(H,1)~=size(H,2))`.

## Parameters / inputs

- parameters.pulse_dur -duration of each pulse, a two-
- element vector, seconds
- parameters.pulse_amp -amplitude of each pulse, a two-
- element vector, rad/s
- parameters.mq_order -MQMAS coherence order
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+
- + the parameters required by the singlerot.m
- context function that will provide H, R, and K

## Outputs

- fid -2D amplitude mode free induction decay

## Header notes

- parameters.sweep should be a two-element vector
- with both elements set equal to parameters.rate,
- this is because this pulse sequence is strobo-
- scopic with respect to the rotor period

## Implementation structure

- Rotor-synchronous MQMAS pulse sequence. Syntax:
- fid=mqmas(spin_system,parameters,H,R,K)
- This function should normally be called using singlerot.m context
- that would provide H, R, and K.
- parameters.pulse_dur -duration of each pulse, a two-
- element vector, seconds
- parameters.pulse_amp -amplitude of each pulse, a two-
- element vector, rad/s
- parameters.mq_order -MQMAS coherence order
- parameters.rho0 -initial condition, usually Lz
- parameters.coil -detection state, usually L+
- + the parameters required by the singlerot.m

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `decouple()`, `step()`, `coherence()`, `evolution()`, `ismatrix()`, `isequal()`, `isfield()`, `iscell()`, `all()`, `cellfun()`, `ismember()`, `any()`, `isscalar()`.
