# experiments/nmr_liquids/roesy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/roesy.m`
- Signature: `fid=roesy(spin_system,parameters,H,R,K)`
- Total lines: 142

## Purpose

Phase-sensitive homonuclear ROESY pulse sequence, assuming ideal spin-lock, described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 47-48: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 50-51: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 53-54: Coherent evolution timestep; implemented by `timestep=1./parameters.sweep`.
- Lines 56-57: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 59-60: Pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 63-64: First pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 66-67: F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 69-70: Analytical spin lock; implemented by `rho_stack_cos=spinlock(spin_system,Lx,Ly,rho_stack,'Y')`.
- Lines 73-74: Mixing time; implemented by `rho_stack_cos=evolution(spin_system,1i*R+1i*K,[],rho_stack_cos,parameters.tmix,1,'final')`.
- Lines 77-78: F2 evolution; implemented by `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 51: computes `L` using `L=H+1i*R+1i*K`.
- Lines 54: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 57: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 60: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 61: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 64: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 67: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 70: computes `rho_stack_cos` using `rho_stack_cos=spinlock(spin_system,Lx,Ly,rho_stack,'Y')`.
- Lines 71: computes `rho_stack_sin` using `rho_stack_sin=spinlock(spin_system,Lx,Ly,rho_stack,'X')`.
- Lines 78: computes `fid.cos` using `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.
- Lines 79: computes `fid.sin` using `fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 84: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Syntax

```matlab
fid=roesy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep -a vector with sweep widths
- in F1 and F2 directions, Hz
- parameters.npoints -a vector with point count
- in F1 and F2 directions
- parameters.spins -nuclei on which the sequence
- runs, e.g. {'1H'}
- parameters.tmix -mixing time, seconds
- parameters.rho0 -initial state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos, fid.sin -components of the free induction
- decay for hypercomplex processing
- Note: this ideal spin-lock model does not represent finite RF
- amplitude, RF offset, Hartmann-Hahn matching errors, or
- explicit RF phase transients.

## Implementation structure

- Phase-sensitive homonuclear ROESY pulse sequence, assuming ideal
- spin-lock, described in:
- fid=roesy(spin_system,parameters,H,R,K)
- parameters.sweep -a vector with sweep widths
- in F1 and F2 directions, Hz
- parameters.npoints -a vector with point count
- in F1 and F2 directions
- parameters.spins -nuclei on which the sequence
- runs, e.g. {'1H'}
- parameters.tmix -mixing time, seconds
- parameters.rho0 -initial state
- H -Hamiltonian matrix, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `timestep()`, `spinlock()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
