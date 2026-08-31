# experiments/nmr_liquids/ct_hsqc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/ct_hsqc.m`
- Signature: `fid=ct_hsqc(spin_system,parameters,H,R,K)`
- Total lines: 222

## Purpose

Constant-time phase-sensitive HSQC pulse sequence from:

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

- Lines 44-45: Set default detection decoupling; implemented by `if ~isfield(parameters,'decouple_f2')`.
- Lines 49-50: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 52-53: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 55-56: Coherent evolution timesteps; implemented by `timestep=1./parameters.sweep`.
- Lines 58-59: J-coupling evolution time; implemented by `tau=abs(1/(4*parameters.J))`.
- Lines 61-62: Initial condition; implemented by `if ~isfield(parameters,'rho0')`.
- Lines 66-67: Detection state; implemented by `if ~isfield(parameters,'coil')`.
- Lines 71-72: Pulse operators; implemented by `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 76-77: Pulse on I spin (1H); implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 79-80: tau evolution; implemented by `rho=evolution(spin_system,L,[],rho,tau,1,'final')`.
- Lines 82-83: Inversion pulses; implemented by `rho=step(spin_system,Hx+Cx,rho,pi)`.
- Lines 88-89: Y pulse on I spin; implemented by `rho=step(spin_system,Hy,rho,pi/2)`.
- Lines 91-92: Pulse on S spin (13C) with coherence selection; implemented by `rho=step(spin_system,Cx,rho,pi/2)-step(spin_system,Cx,rho,-pi/2)`.
- Lines 94-95: Preallocate rho stack; implemented by `rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i)`.
- Lines 97-98: Get the time grid for the CT period; implemented by `t1_grid=(0:(parameters.npoints(1)-1))/parameters.sweep(1)`.
- Lines 101-102: Loop over the value of t1_grid; implemented by `for n=1:parameters.npoints(1)`.
- Lines 104-105: Run the delay; implemented by `rho_current=evolution(spin_system,L,[],rho,(CT-t1_grid(n))/2,1,'final')`.
- Lines 107-108: Run the pulse on the S spin; implemented by `rho_current=step(spin_system,Cx,rho_current,pi)`.

### Control flow inferred from the code

- Line 45: conditional branch on `~isfield(parameters,'decouple_f2')`.
- Line 62: conditional branch on `~isfield(parameters,'rho0')`.
- Line 67: conditional branch on `~isfield(parameters,'coil')`.
- Line 102: `for` loop over `n=1:parameters.npoints(1)`.

### Key state/data transformations

- Lines 46: computes `parameters.decouple_f2` using `parameters.decouple_f2={}`.
- Lines 53: computes `L` using `L=H+1i*R+1i*K`.
- Lines 56: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 59: computes `tau` using `tau=abs(1/(4*parameters.J))`.
- Lines 63: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 68: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 72: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 73: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 74: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 77: computes `rho` using `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 95: computes `rho_stack` using `rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i)`.
- Lines 98: computes `t1_grid` using `t1_grid=(0:(parameters.npoints(1)-1))/parameters.sweep(1)`.
- Lines 99: computes `CT` using `CT=t1_grid(end)`.
- Lines 105: computes `rho_current` using `rho_current=evolution(spin_system,L,[],rho,(CT-t1_grid(n))/2,1,'final')`.
- Lines 120: computes `rho_stack(:,n)` using `rho_stack(:,n)=rho_current`.
- Lines 125-126: computes `rho_stack_pos` using `rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0}, {parameters.spins{1},+1}})`.
- Lines 127-128: computes `rho_stack_neg` using `rho_stack_neg=coherence(spin_system,rho_stack,{{parameters.spins{2}, 0}, {parameters.spins{1},-1}})`.
- Lines 147: computes `[L,rho_stack_pos]` using `[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2)`.

### Local helper functions

- Line 159: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=ct_hsqc(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 [optional] nuclei to decouple
- in F2, e.g. {'15N','13C'}
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

- Constant-time phase-sensitive HSQC pulse sequence from:
- fid=ct_hsqc(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 [optional] nuclei to decouple
- in F2, e.g. {'15N','13C'}
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.pos,fid.neg - two components of the States quadrature

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `t1_grid()`, `rho_stack()`, `coherence()`, `decouple()`, `timestep()`, `ismember()`, `ismatrix()`, `all()`, `elseif()`, `any()`.
