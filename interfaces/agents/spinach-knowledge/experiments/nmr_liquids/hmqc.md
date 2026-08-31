# experiments/nmr_liquids/hmqc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/hmqc.m`
- Signature: `fid=hmqc(spin_system,parameters,H,R,K)`
- Total lines: 184

## Purpose

Magnitude-mode HMQC pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 49-50: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 52-53: Evolution timesteps; implemented by `timestep=1./parameters.sweep`.
- Lines 55-56: J-coupling evolution time; implemented by `delta=abs(1/(2*parameters.J))`.
- Lines 58-59: Initial state; implemented by `rho=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 61-62: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 64-65: Pulse operators; implemented by `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 68-69: Pulse on F2; implemented by `rho=step(spin_system,Lx_F2,rho,pi/2)`.
- Lines 71-72: Delta evolution; implemented by `rho=evolution(spin_system,L,[],rho,delta,1,'final')`.
- Lines 74-75: Pulse on F1; implemented by `rho=step(spin_system,Lx_F1,rho,pi/2)`.
- Lines 77-78: Coherence selection; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 80-82: First half of F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 84-85: F1 dimension decoupling; implemented by `for n=1:numel(parameters.decouple_f1)`.
- Lines 90-92: Second half of F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2, parameters.npoints(1)-1,'refocus')`.
- Lines 94-95: Pulse on F1; implemented by `rho_stack=step(spin_system,Lx_F1,rho_stack,pi/2)`.
- Lines 97-98: Delta evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,delta,1,'final')`.
- Lines 100-101: F2 dimension decoupling; implemented by `[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.decouple_f2)`.
- Lines 103-105: Detection on F2; implemented by `fid=evolution(spin_system,L,coil,rho_stack,timestep(2), parameters.npoints(2)-1,'observable')`.

### Control flow inferred from the code

- Line 85: `for` loop over `n=1:numel(parameters.decouple_f1)`.

### Key state/data transformations

- Lines 50: computes `L` using `L=H+1i*R+1i*K`.
- Lines 53: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 56: computes `delta` using `delta=abs(1/(2*parameters.J))`.
- Lines 59: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 62: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 65: computes `Lx_F1` using `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 66: computes `Lx_F2` using `Lx_F2=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 81-82: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 86: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.decouple_f1{n})`.
- Lines 101: computes `[L,rho_stack]` using `[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.decouple_f2)`.
- Lines 104-105: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,timestep(2), parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 110: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Syntax

```matlab
fid=hmqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths in each
- frequency direction, Hz
- parameters.npoints [F1 F2] numbers of points in
- each time direction
- parameters.spins {F1 F2} nuclei, e.g. {'15N','1H'}
- parameters.decouple_f2 nuclei to decouple in F2,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint 180-degree
- refocusing pulses in F1, e.g. {'1H'}
- parameters.J primary scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay for magnitude mode processing
- Note: natural abundance experiments should make use of the iso-
- tope dilution functionality. See dilute.m function.

## Implementation structure

- Magnitude-mode HMQC pulse sequence from:
- fid=hmqc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths in each
- frequency direction, Hz
- parameters.npoints [F1 F2] numbers of points in
- each time direction
- parameters.spins {F1 F2} nuclei, e.g. {'15N','1H'}
- parameters.decouple_f2 nuclei to decouple in F2,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint 180-degree
- refocusing pulses in F1, e.g. {'1H'}
- parameters.J primary scalar coupling, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `timestep()`, `decouple()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
