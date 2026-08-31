# experiments/nmr_liquids/inadequate_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/inadequate_2d.m`
- Signature: `fid=inadequate_2d(spin_system,parameters,H,R,K)`
- Total lines: 166

## Purpose

2D INADEQUATE pulse sequence from:

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

- Lines 42-43: Set default decoupling; implemented by `if ~isfield(parameters,'decouple')`.
- Lines 47-48: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 50-51: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 53-54: Decoupling; implemented by `L=decouple(spin_system,L,[],parameters.decouple)`.
- Lines 56-57: Get timing parameters; implemented by `tau=abs(1/(4*parameters.J))`.
- Lines 60-61: Initial state along Z; implemented by `rho=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 63-64: Quadrature detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1})`.
- Lines 66-67: Pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 71-72: 90x pulse; implemented by `rho_cos=step(spin_system,Lx,rho,pi/2)`.
- Lines 75-76: J-coupling evolution during tau; implemented by `rho_cos=evolution(spin_system,L,[],rho_cos,tau,1,'final')`.
- Lines 79-80: Inversion pulse; implemented by `rho_cos=step(spin_system,Lx,rho_cos,pi)`.
- Lines 83-84: Second tau evolution; implemented by `rho_cos=evolution(spin_system,L,[],rho_cos,tau,1,'final')`.
- Lines 87-88: 90x pulse, creation of MQC; implemented by `rho_cos=step(spin_system,Lx,rho_cos,pi/2)`.
- Lines 91-92: Apply the double-quantum filter; implemented by `rho_cos=coherence(spin_system,rho_cos,{{parameters.spins{1},[+2,-2]}})`.
- Lines 95-96: t1 evolution; implemented by `rho_stack_cos=evolution(spin_system,L,[],rho_cos,timestep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 99-100: 90x pulse, creation of observable magnetization; implemented by `rho_stack_cos=step(spin_system,Lx,rho_stack_cos,pi/2)`.
- Lines 103-104: Run the t2 evolution; implemented by `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.

### Control flow inferred from the code

- Line 43: conditional branch on `~isfield(parameters,'decouple')`.

### Key state/data transformations

- Lines 44: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 51: computes `L` using `L=H+1i*R+1i*K`.
- Lines 57: computes `tau` using `tau=abs(1/(4*parameters.J))`.
- Lines 58: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 61: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 64: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1})`.
- Lines 67: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 68: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 69: computes `Lxy` using `Lxy=cos(pi/4)*Lx+sin(pi/4)*Ly`.
- Lines 72: computes `rho_cos` using `rho_cos=step(spin_system,Lx,rho,pi/2)`.
- Lines 73: computes `rho_sin` using `rho_sin=step(spin_system,Lxy,rho,pi/2)`.
- Lines 96: computes `rho_stack_cos` using `rho_stack_cos=evolution(spin_system,L,[],rho_cos,timestep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 97: computes `rho_stack_sin` using `rho_stack_sin=evolution(spin_system,L,[],rho_sin,timestep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 104: computes `fid.cos` using `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.
- Lines 105: computes `fid.sin` using `fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 110: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=inadequate_2d(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins active nucleus, e.g. {'13C'}
- parameters.decouple optional nuclei to decouple, e.g. {'1H'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos,fid.sin -two components of the States signal
- Notes: use dilute.m to generate carbon pair isotopomers. The F1
- axis is a double-quantum frequency coordinate.
- Theresa Hune
- Christian Griesinger

## Implementation structure

- 2D INADEQUATE pulse sequence from:
- fid=inadequate_2d(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins active nucleus, e.g. {'13C'}
- parameters.decouple optional nuclei to decouple, e.g. {'1H'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.cos,fid.sin - two components of the States signal
- axis is a double-quantum frequency coordinate.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `decouple()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `timestep()`, `ismember()`, `ismatrix()`, `all()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
