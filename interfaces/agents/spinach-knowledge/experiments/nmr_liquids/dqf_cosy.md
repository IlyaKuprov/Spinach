# experiments/nmr_liquids/dqf_cosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/dqf_cosy.m`
- Signature: `fid=dqf_cosy(spin_system,parameters,H,R,K)`
- Total lines: 128

## Purpose

Phase-sensitive double-quantum filtered COSY pulse sequence. Implemented as described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 42-43: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 45-46: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 48-49: Compute the evolution timestep; implemented by `timestep=1/parameters.sweep`.
- Lines 51-52: Initial state; implemented by `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 54-55: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 57-58: Get the pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 61-62: Apply the first pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 64-65: Run the F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints(1)-1,'trajectory')`.
- Lines 67-68: Apply the second pulse (States hypercomplex); implemented by `rho_stack_sin=step(spin_system,Lx,rho_stack,pi/2)`.
- Lines 71-72: Apply the double-quantum filter; implemented by `rho_stack_cos=coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+2,-2]}})`.
- Lines 75-76: Apply the third pulse; implemented by `rho_stack_sin=step(spin_system,Lx,rho_stack_sin,pi/2)`.
- Lines 79-80: Run the F2 evolution; implemented by `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep,parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 46: computes `L` using `L=H+1i*R+1i*K`.
- Lines 49: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 52: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 55: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 58: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 59: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 65: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep,parameters.npoints(1)-1,'trajectory')`.
- Lines 68: computes `rho_stack_sin` using `rho_stack_sin=step(spin_system,Lx,rho_stack,pi/2)`.
- Lines 69: computes `rho_stack_cos` using `rho_stack_cos=step(spin_system,Ly,rho_stack,pi/2)`.
- Lines 80: computes `fid.cos` using `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep,parameters.npoints(2)-1,'observable')`.
- Lines 81: computes `fid.sin` using `fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep,parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 86: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=dqf_cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos, fid.sin -components of the free induction
- decay for hypercomplex processing
- Note: the double-quantum filter is implemented as an exact
- analytical coherence-order projection, not as an explicit
- phase cycle or finite-gradient selection block.

## Implementation structure

- Phase-sensitive double-quantum filtered COSY pulse sequence.
- Implemented as described in:
- fid=dqf_cosy(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.cos, fid.sin - components of the free induction
- decay for hypercomplex processing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `nnz()`, `strcmp()`, `any()`.
