# experiments/nmr_liquids/pansy_triple.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/pansy_triple.m`
- Signature: `fid=pansy_triple(spin_system,parameters,H,R,K)`
- Total lines: 149

## Purpose

Triple-channel PANSY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 49-50: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 52-53: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 55-56: Timestep; implemented by `timesteps=1./parameters.sweep`.
- Lines 58-59: Initial state; implemented by `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 61-62: Detection state; implemented by `coil_h=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 66-67: Get pulse operators; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 72-73: Apply the first pulse; implemented by `rho_p=step(spin_system,Hx,rho,+pi/2)`.
- Lines 76-77: Select "+1" coherence; implemented by `rho_p=coherence(spin_system,rho_p,{{parameters.spins{1},+1}})`.
- Lines 80-82: Run F1 evolution; implemented by `rho_stack_p=evolution(spin_system,L,[],rho_p,timesteps(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 86-87: Apply the second pulse pair; implemented by `rho_stack_p=step(spin_system,Hy+X1y+X2y,rho_stack_p,pi/2)`.
- Lines 90-91: Perform axial peak elimination; implemented by `rho_stack=rho_stack_p-rho_stack_m`.
- Lines 93-95: Run F2 evolution; implemented by `fid.aa=evolution(spin_system,L,coil_h,rho_stack,timesteps(1), parameters.npoints(1)-1,'observable')`.

### Key state/data transformations

- Lines 53: computes `L` using `L=H+1i*R+1i*K`.
- Lines 56: computes `timesteps` using `timesteps=1./parameters.sweep`.
- Lines 59: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 62: computes `coil_h` using `coil_h=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 63: computes `coil_x1` using `coil_x1=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 64: computes `coil_x2` using `coil_x2=state(spin_system,'L+',parameters.spins{3},'cheap')`.
- Lines 67: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 68: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 69: computes `X1y` using `X1y=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 70: computes `X2y` using `X2y=operator(spin_system,'Ly',parameters.spins{3})`.
- Lines 73: computes `rho_p` using `rho_p=step(spin_system,Hx,rho,+pi/2)`.
- Lines 74: computes `rho_m` using `rho_m=step(spin_system,Hx,rho,-pi/2)`.
- Lines 81-82: computes `rho_stack_p` using `rho_stack_p=evolution(spin_system,L,[],rho_p,timesteps(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 83-84: computes `rho_stack_m` using `rho_stack_m=evolution(spin_system,L,[],rho_m,timesteps(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 91: computes `rho_stack` using `rho_stack=rho_stack_p-rho_stack_m`.
- Lines 94-95: computes `fid.aa` using `fid.aa=evolution(spin_system,L,coil_h,rho_stack,timesteps(1), parameters.npoints(1)-1,'observable')`.
- Lines 96-97: computes `fid.ab` using `fid.ab=evolution(spin_system,L,coil_x1,rho_stack,timesteps(2), parameters.npoints(2)-1,'observable')`.
- Lines 98-99: computes `fid.ac` using `fid.ac=evolution(spin_system,L,coil_x2,rho_stack,timesteps(3), parameters.npoints(3)-1,'observable')`.

### Local helper functions

- Line 104: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=pansy_triple(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -three nuclei on which the sequence
- runs, e.g. {'1H','13C','15N'}
- parameters.sweep -a vector with three sweep widths
- in Hz
- parameters.npoints -a vector of three integers specify-
- ing point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.aa -COSY FID on F1 nuclei
- fid.ab -COSY FID on F1,F2 nuclei
- fid.ac -COSY FID on F1,F3 nuclei
- Note: decoupling with respect to any of the working nuclei is
- impossible in this pulse sequence.
- Note: this is the magnitude-mode analytical-pathway version.
- Gradient echo/anti-echo PANSY variants should be implemented
- as separate pulse sequence functions.
- Andrew Porter

## Implementation structure

- Triple-channel PANSY pulse sequence from:
- fid=pansy_triple(spin_system,parameters,H,R,K)
- parameters.spins -three nuclei on which the sequence
- runs, e.g. {'1H','13C','15N'}
- parameters.sweep -a vector with three sweep widths
- in Hz
- parameters.npoints -a vector of three integers specify-
- ing point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.aa -COSY FID on F1 nuclei

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `coherence()`, `evolution()`, `timesteps()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
