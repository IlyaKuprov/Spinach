# experiments/nmr_liquids/hmbc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/hmbc.m`
- Signature: `fid=hmbc(spin_system,parameters,H,R,K)`
- Total lines: 168

## Purpose

Magnitude-mode HMBC pulse sequence from:

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
- Lines 52-53: Get evolution time step; implemented by `timestep=1./parameters.sweep`.
- Lines 55-56: Fixed evolution times; implemented by `delta_a=abs(1/(2*parameters.J))`.
- Lines 59-60: Initial and detection states; implemented by `rho0=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 63-64: Pulse operators; implemented by `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 67-68: First proton pulse; implemented by `rho=step(spin_system,Hx,rho0,pi/2)`.
- Lines 70-71: First evolution period; implemented by `rho=evolution(spin_system,L,[],rho,delta_a,1,'final')`.
- Lines 73-74: First carbon pulse; implemented by `rho=step(spin_system,Cx,rho,+pi/2)`.
- Lines 76-77: Second evolution period; implemented by `rho=evolution(spin_system,L,[],rho,delta_b,1,'final')`.
- Lines 79-81: Second carbon pulse; implemented by `rho=step(spin_system,Cx,rho,+pi/2)- step(spin_system,Cx,rho,-pi/2)`.
- Lines 83-84: Coherence selection; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 86-88: F1 evolution period, first half; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 90-91: Proton decoupling pulse; implemented by `rho_stack=step(spin_system,Hx,rho_stack,pi)`.
- Lines 93-95: F1 evolution period, second half; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2, parameters.npoints(1)-1,'refocus')`.
- Lines 97-98: Third carbon pulse; implemented by `rho_stack=step(spin_system,Cx,rho_stack,pi/2)`.
- Lines 100-102: Detection; implemented by `fid=evolution(spin_system,L,coil,rho_stack,timestep(2), parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 50: computes `L` using `L=H+1i*R+1i*K`.
- Lines 53: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 56: computes `delta_a` using `delta_a=abs(1/(2*parameters.J))`.
- Lines 57: computes `delta_b` using `delta_b=parameters.delta_b`.
- Lines 60: computes `rho0` using `rho0=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 61: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 64: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 65: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 68: computes `rho` using `rho=step(spin_system,Hx,rho0,pi/2)`.
- Lines 87-88: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 101-102: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,timestep(2), parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 107: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Syntax

```matlab
fid=hmbc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths in each
- frequency direction, Hz
- parameters.npoints [F1 F2] numbers of points in
- each time direction
- parameters.spins {F1 F2} nuclei, e.g. {'15N','1H'}
- parameters.J primary scalar coupling, Hz
- parameters.delta_b delta_2 delay from the paper
- cited above; the authors recom-
- mend 60e-3 seconds
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay for magnitude mode processing
- Note: natural abundance experiments should make use of the iso-
- tope dilution functionality. See dilute.m function.

## Implementation structure

- Magnitude-mode HMBC pulse sequence from:
- fid=hmbc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths in each
- frequency direction, Hz
- parameters.npoints [F1 F2] numbers of points in
- each time direction
- parameters.spins {F1 F2} nuclei, e.g. {'15N','1H'}
- parameters.J primary scalar coupling, Hz
- parameters.delta_b delta_2 delay from the paper
- cited above; the authors recom-
- mend 60e-3 seconds
- H -Hamiltonian matrix, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `timestep()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`, `strcmp()`.
