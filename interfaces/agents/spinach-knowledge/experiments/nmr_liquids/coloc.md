# experiments/nmr_liquids/coloc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/coloc.m`
- Signature: `fid=coloc(spin_system,parameters,H,R,K)`
- Total lines: 187

## Purpose

COLOC NMR pulse sequence from: Implemented as shown in Fig 1b, without the dashed pulses during the delta(2) period. Delta(1) defaults to half of the maximum F1 evolution time implied by the sweep width, delta(2) must be spe- cified. Syntax: fid=coloc(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 49-50: Set default delta1 delay; implemented by `if ~isfield(parameters,'delta1')`.
- Lines 54-55: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 57-58: Sequence timing; implemented by `timestep=1./parameters.sweep`.
- Lines 61-62: Initial state; implemented by `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 64-65: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 67-68: Pulse operators; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 72-73: 90 degree pulse on H; implemented by `rho=step(spin_system,Hx,rho,pi/2)`.
- Lines 75-77: Report the delta(1); implemented by `report(spin_system,['COLOC delta(1), seconds: ' num2str(parameters.delta1)])`.
- Lines 79-80: Preallocate the evolution stack; implemented by `rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i)`.
- Lines 82-83: Loop over the embedded-echo evolution time; implemented by `for n=1:parameters.npoints(1)`.
- Lines 85-87: First part of the echo block; implemented by `rho_current=evolution(spin_system,L,[],rho, t1_grid(n)/2,1,'final')`.
- Lines 89-90: Echo pulse pair; implemented by `rho_current=step(spin_system,Hx,rho_current,pi)`.
- Lines 93-95: Second part of the echo block; implemented by `rho_current=evolution(spin_system,L,[],rho_current, parameters.delta1-t1_grid(n)/2,1,'final')`.
- Lines 97-98: Assign the stack element; implemented by `rho_stack(:,n)=rho_current`.
- Lines 102-103: Proton coherence selection; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},-1}})`.
- Lines 105-106: 90 degree pulses on H and C; implemented by `rho_stack=step(spin_system,Hx,rho_stack,pi/2)`.
- Lines 109-110: Delta evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,parameters.delta2,1,'final')`.

### Control flow inferred from the code

- Line 50: conditional branch on `~isfield(parameters,'delta1')`.
- Line 83: `for` loop over `n=1:parameters.npoints(1)`.

### Key state/data transformations

- Lines 51: computes `parameters.delta1` using `parameters.delta1=(parameters.npoints(1)-1)/(2*parameters.sweep(1))`.
- Lines 55: computes `L` using `L=H+1i*R+1i*K`.
- Lines 58: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 59: computes `t1_grid` using `t1_grid=(0:(parameters.npoints(1)-1))*timestep(1)`.
- Lines 62: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 65: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 68: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 69: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 70: computes `Cy` using `Cy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 80: computes `rho_stack` using `rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i)`.
- Lines 86-87: computes `rho_current` using `rho_current=evolution(spin_system,L,[],rho, t1_grid(n)/2,1,'final')`.
- Lines 98: computes `rho_stack(:,n)` using `rho_stack(:,n)=rho_current`.
- Lines 116: computes `[L,rho_stack]` using `[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.spins(1))`.
- Lines 119-120: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,timestep(2), parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 125: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.delta2 COLOC delta2 (see the paper),
- typically 40e-3 seconds
- parameters.delta1 optional COLOC delta1 delay, seconds;
- must be at least half of maximum t1
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay for magnitude mode processing
- Note: natural abundance simulations should make use of the isotope
- dilution functionality. See dilute.m function.

## Implementation structure

- COLOC NMR pulse sequence from:
- Implemented as shown in Fig 1b, without the dashed pulses during
- the delta(2) period. Delta(1) defaults to half of the maximum F1
- evolution time implied by the sweep width, delta(2) must be spe-
- cified. Syntax:
- fid=coloc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.delta2 COLOC delta2 (see the paper),
- typically 40e-3 seconds
- parameters.delta1 optional COLOC delta1 delay, seconds;

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `timestep()`, `state()`, `operator()`, `step()`, `report()`, `delta()`, `num2str()`, `evolution()`, `t1_grid()`, `rho_stack()`, `coherence()`, `decouple()`, `ismember()`, `ismatrix()`.
