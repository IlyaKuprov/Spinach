# experiments/nmr_liquids/hsqc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/hsqc.m`
- Signature: `fid=hsqc(spin_system,parameters,H,R,K)`
- Total lines: 222

## Purpose

Phase-sensitive HSQC pulse sequence from:

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

- Lines 48-49: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 51-52: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 54-55: Coherent evolution timesteps; implemented by `timestep=1./parameters.sweep`.
- Lines 57-58: J-coupling evolution time; implemented by `delta=abs(1/(2*parameters.J))`.
- Lines 60-61: Initial condition; implemented by `if ~isfield(parameters,'rho0')`.
- Lines 65-66: Detection state; implemented by `if ~isfield(parameters,'coil')`.
- Lines 70-71: Pulse operators; implemented by `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 75-76: Pulse on F2; implemented by `rho=step(spin_system,Lx_F2,parameters.rho0,pi/2)`.
- Lines 78-79: Delta evolution; implemented by `rho=evolution(spin_system,L,[],rho,delta/2,1,'final')`.
- Lines 81-82: Inversion pulses; implemented by `rho=step(spin_system,Lx_F2+Lx_F1,rho,pi)`.
- Lines 87-88: Pulse on F2; implemented by `rho=step(spin_system,Ly_F2,rho,pi/2)`.
- Lines 90-91: Pulse on F1 with coherence selection; implemented by `rho=step(spin_system,Lx_F1,rho,pi/2)-step(spin_system,Lx_F1,rho,-pi/2)`.
- Lines 93-95: F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 103-105: Coherence selection; implemented by `rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0}, {parameters.spins{1},+1}})`.
- Lines 109-110: Pulses on F1 and F2; implemented by `rho_stack_pos=step(spin_system,Lx_F1+Lx_F2,rho_stack_pos,pi/2)`.
- Lines 113-114: Delta evolution; implemented by `rho_stack_pos=evolution(spin_system,L,[],rho_stack_pos,delta/2,1,'final')`.
- Lines 117-118: Pulses on F1 and F2; implemented by `rho_stack_pos=step(spin_system,Lx_F1+Lx_F2,rho_stack_pos,pi)`.
- Lines 125-126: Decoupling in F2; implemented by `[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2)`.

### Control flow inferred from the code

- Line 61: conditional branch on `~isfield(parameters,'rho0')`.
- Line 66: conditional branch on `~isfield(parameters,'coil')`.
- Line 96: `for` loop over `n=1:numel(parameters.decouple_f1)`.

### Key state/data transformations

- Lines 52: computes `L` using `L=H+1i*R+1i*K`.
- Lines 55: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 58: computes `delta` using `delta=abs(1/(2*parameters.J))`.
- Lines 62: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 67: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 71: computes `Lx_F1` using `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 72: computes `Lx_F2` using `Lx_F2=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 73: computes `Ly_F2` using `Ly_F2=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 76: computes `rho` using `rho=step(spin_system,Lx_F2,parameters.rho0,pi/2)`.
- Lines 94-95: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 97: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.decouple_f1{n})`.
- Lines 104-105: computes `rho_stack_pos` using `rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0}, {parameters.spins{1},+1}})`.
- Lines 106-107: computes `rho_stack_neg` using `rho_stack_neg=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0}, {parameters.spins{1},-1}})`.
- Lines 126: computes `[L,rho_stack_pos]` using `[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2)`.
- Lines 127: computes `[L,rho_stack_neg]` using `[L,rho_stack_neg]=decouple(spin_system,L,rho_stack_neg,parameters.decouple_f2)`.
- Lines 130-131: computes `fid.pos` using `fid.pos=evolution(spin_system,L,parameters.coil,rho_stack_pos, timestep(2),parameters.npoints(2)-1,'observable')`.
- Lines 132-133: computes `fid.neg` using `fid.neg=evolution(spin_system,L,parameters.coil,rho_stack_neg, timestep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 138: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=hsqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 nuclei to decouple in F2, e.g.
- {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint
- 180-degree refocusing pulses in
- F1, e.g. {'1H','13C'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.pos,fid.neg -two components of the States quadrature
- signal.
- Note: natural abundance simulations should make use of the isotope
- dilution functionality. See dilute.m function.

## Implementation structure

- Phase-sensitive HSQC pulse sequence from:
- fid=hsqc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 nuclei to decouple in F2, e.g.
- {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint
- 180-degree refocusing pulses in
- F1, e.g. {'1H','13C'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `state()`, `operator()`, `step()`, `evolution()`, `timestep()`, `coherence()`, `decouple()`, `ismember()`, `ismatrix()`, `all()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
