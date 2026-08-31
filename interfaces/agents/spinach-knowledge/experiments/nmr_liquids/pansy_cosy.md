# experiments/nmr_liquids/pansy_cosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/pansy_cosy.m`
- Signature: `fid=pansy_cosy(spin_system,parameters,H,R,K)`
- Total lines: 144

## Purpose

Magnitude mode PANSY-COSY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 49-50: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 52-53: Get evolution time steps; implemented by `timesteps=1./parameters.sweep`.
- Lines 55-56: Initial state; implemented by `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 58-59: Detection state; implemented by `coil_h=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 62-63: Get pulse operators; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 67-68: Apply the first pulse; implemented by `rho_p=step(spin_system,Hx,rho,+pi/2)`.
- Lines 71-72: Select "+1" coherence; implemented by `rho_p=coherence(spin_system,rho_p,{{parameters.spins{1},+1}})`.
- Lines 75-77: Run F1 evolution; implemented by `rho_stack_p=evolution(spin_system,L,[],rho_p,timesteps(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 81-82: Apply the second pulse pair; implemented by `rho_stack_p=step(spin_system,Hy+Cy,rho_stack_p,pi/2)`.
- Lines 85-86: Perform axial peak elimination; implemented by `rho_stack=rho_stack_p-rho_stack_m`.
- Lines 88-90: Run F2 evolution; implemented by `fid.aa=evolution(spin_system,L,coil_h,rho_stack,timesteps(1), parameters.npoints(1)-1,'observable')`.

### Key state/data transformations

- Lines 50: computes `L` using `L=H+1i*R+1i*K`.
- Lines 53: computes `timesteps` using `timesteps=1./parameters.sweep`.
- Lines 56: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 59: computes `coil_h` using `coil_h=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 60: computes `coil_c` using `coil_c=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 63: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 64: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 65: computes `Cy` using `Cy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 68: computes `rho_p` using `rho_p=step(spin_system,Hx,rho,+pi/2)`.
- Lines 69: computes `rho_m` using `rho_m=step(spin_system,Hx,rho,-pi/2)`.
- Lines 76-77: computes `rho_stack_p` using `rho_stack_p=evolution(spin_system,L,[],rho_p,timesteps(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 78-79: computes `rho_stack_m` using `rho_stack_m=evolution(spin_system,L,[],rho_m,timesteps(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 86: computes `rho_stack` using `rho_stack=rho_stack_p-rho_stack_m`.
- Lines 89-90: computes `fid.aa` using `fid.aa=evolution(spin_system,L,coil_h,rho_stack,timesteps(1), parameters.npoints(1)-1,'observable')`.
- Lines 91-92: computes `fid.ab` using `fid.ab=evolution(spin_system,L,coil_c,rho_stack,timesteps(2), parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 97: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=pansy_cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.spins -nuclei on which the sequence runs,
- e.g. {'1H','13C'}
- parameters.sweep -a vector with two sweep widths in Hz
- parameters.npoints -a vector of integers specifying
- point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.aa -magnitude mode COSY FID on F1,F1 nuclei
- fid.ab -magnitude mode COSY FID on F1,F2 nuclei
- Note: decoupling with respect to either of the working nuclei
- is impossible in this pulse sequence.
- Note: this is the magnitude-mode analytical-pathway version.
- Gradient echo/anti-echo PANSY-COSY variants should be
- implemented as separate pulse sequence functions.
- Andrew Porter

## Implementation structure

- Magnitude mode PANSY-COSY pulse sequence from:
- fid=pansy_cosy(spin_system,parameters,H,R,K)
- parameters.spins -nuclei on which the sequence runs,
- e.g. {'1H','13C'}
- parameters.sweep -a vector with two sweep widths in Hz
- parameters.npoints -a vector of integers specifying
- point count in each dimension
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.aa -magnitude mode COSY FID on F1,F1 nuclei
- fid.ab -magnitude mode COSY FID on F1,F2 nuclei

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `coherence()`, `evolution()`, `timesteps()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
