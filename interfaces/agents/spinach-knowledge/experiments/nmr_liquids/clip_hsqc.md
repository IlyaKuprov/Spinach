# experiments/nmr_liquids/clip_hsqc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/clip_hsqc.m`
- Signature: `fid=clip_hsqc(spin_system,parameters,H,R,K)`
- Total lines: 188

## Purpose

CLIP-HSQC pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 41-42: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 44-45: Coherent evolution timesteps; implemented by `timestep=1./parameters.sweep`.
- Lines 47-48: J-coupling evolution time; implemented by `delta=abs(1/(2*parameters.J))`.
- Lines 50-51: Initial state; implemented by `rho=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 53-54: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 56-57: Pulse operators; implemented by `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 61-62: First pulse; implemented by `report(spin_system,'propagating initial state forwards in time ')`.
- Lines 65-66: Delta evolution; implemented by `rho=evolution(spin_system,L,[],rho,delta/2,1,'final')`.
- Lines 68-69: Second pulse; implemented by `rho=step(spin_system,Lx_F2+Lx_F1,rho,pi)`.
- Lines 74-75: Pulses; implemented by `rho=step(spin_system,Ly_F2,rho,pi/2)`.
- Lines 79-80: F1 evolution; implemented by `rho_stack_px=evolution(spin_system,L,[],rho_px,timestep(1)/2,parameters.npoints(1)-1,'trajectory')`.
- Lines 87-88: Coherence selection by the first gradient; implemented by `rho_stack_px_P=coherence(spin_system,rho_stack_px,{{parameters.spins{2},0},{parameters.spins{1},+1}})`.
- Lines 93-94: F2 evolution; implemented by `report(spin_system,'propagating coil state backwards in time ')`.
- Lines 97-98: Pulses; implemented by `coil_stack_px=step(spin_system,+Lx_F1,coil_stack,-pi/2)`.
- Lines 101-102: Delta evolution; implemented by `coil_stack_px=evolution(spin_system,L',[],coil_stack_px,-delta/2,1,'final')`.
- Lines 105-106: Coherence selection by the second gradient; implemented by `coil_stack_px=coherence(spin_system,coil_stack_px,{{parameters.spins{1},0},{parameters.spins{2},+1}})`.
- Lines 109-110: Pulses; implemented by `coil_stack_px=step(spin_system,Ly_F2+Lx_F1,coil_stack_px,-pi)`.

### Key state/data transformations

- Lines 42: computes `L` using `L=H+1i*R+1i*K`.
- Lines 45: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 48: computes `delta` using `delta=abs(1/(2*parameters.J))`.
- Lines 51: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 54: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 57: computes `Lx_F1` using `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 58: computes `Lx_F2` using `Lx_F2=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 59: computes `Ly_F2` using `Ly_F2=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 76: computes `rho_px` using `rho_px=step(spin_system,+Lx_F1,rho,pi/2)`.
- Lines 77: computes `rho_mx` using `rho_mx=step(spin_system,-Lx_F1,rho,pi/2)`.
- Lines 80: computes `rho_stack_px` using `rho_stack_px=evolution(spin_system,L,[],rho_px,timestep(1)/2,parameters.npoints(1)-1,'trajectory')`.
- Lines 81: computes `rho_stack_mx` using `rho_stack_mx=evolution(spin_system,L,[],rho_mx,timestep(1)/2,parameters.npoints(1)-1,'trajectory')`.
- Lines 88: computes `rho_stack_px_P` using `rho_stack_px_P=coherence(spin_system,rho_stack_px,{{parameters.spins{2},0},{parameters.spins{1},+1}})`.
- Lines 89: computes `rho_stack_mx_P` using `rho_stack_mx_P=coherence(spin_system,rho_stack_mx,{{parameters.spins{2},0},{parameters.spins{1},+1}})`.
- Lines 90: computes `rho_stack_px_N` using `rho_stack_px_N=coherence(spin_system,rho_stack_px,{{parameters.spins{2},0},{parameters.spins{1},-1}})`.
- Lines 91: computes `rho_stack_mx_N` using `rho_stack_mx_N=coherence(spin_system,rho_stack_mx,{{parameters.spins{2},0},{parameters.spins{1},-1}})`.
- Lines 95: computes `coil_stack` using `coil_stack=evolution(spin_system,L',[],coil,-timestep(2),parameters.npoints(2)-1,'trajectory')`.
- Lines 98: computes `coil_stack_px` using `coil_stack_px=step(spin_system,+Lx_F1,coil_stack,-pi/2)`.

### Local helper functions

- Line 134: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=clip_hsqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. {'13C','1H'})
- parameters.J active scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay. States quadrature detection is
- used; the two components of the States signal are re-
- turned in fid.pos and fid.neg

## Implementation structure

- CLIP-HSQC pulse sequence from:
- fid=clip_hsqc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. {'13C','1H'})
- parameters.J active scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid -free induction decay. States quadrature detection is
- used; the two components of the States signal are re-
- turned in fid.pos and fid.neg

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `report()`, `step()`, `evolution()`, `timestep()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
