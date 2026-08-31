# experiments/nmr_liquids/hoesy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/hoesy.m`
- Signature: `fid=hoesy(spin_system,parameters,H,R,K)`
- Total lines: 199

## Purpose

Phase-sensitive heteronuclear NOESY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 53-54: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 56-57: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 59-60: Coherent evolution timestep; implemented by `timestep=1./parameters.sweep`.
- Lines 62-63: Detection state up to a constant multiplier; implemented by `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 65-66: Pulse operators; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 70-71: First pulse on F1; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 73-75: First half of F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 77-78: F1 dimension decoupling; implemented by `for n=1:numel(parameters.decouple_f1)`.
- Lines 83-85: Second half of F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2, parameters.npoints(1)-1,'refocus')`.
- Lines 87-88: Second pulse on F1; implemented by `rho_stack_cos_p=step(spin_system,Hx,rho_stack,+pi/2)`.
- Lines 93-94: Homospoil; implemented by `rho_stack_cos_p=homospoil(spin_system,rho_stack_cos_p,'destroy')`.
- Lines 99-100: Mixing time; implemented by `rho_stack_cos_p=evolution(spin_system,1i*R+1i*K,[],rho_stack_cos_p,parameters.tmix,1,'final')`.
- Lines 111-112: Pulse on F2; implemented by `rho_stack_cos_p=step(spin_system,Cy,rho_stack_cos_p,pi/2)`.
- Lines 117-118: Axial peak elimination in F2; implemented by `rho_stack_cos=rho_stack_cos_p-rho_stack_cos_m`.
- Lines 121-122: Decouple protons; implemented by `[L,rho_stack_cos]=decouple(spin_system,L,rho_stack_cos,parameters.spins(1))`.
- Lines 125-126: F2 evolution; implemented by `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.

### Control flow inferred from the code

- Line 78: `for` loop over `n=1:numel(parameters.decouple_f1)`.

### Key state/data transformations

- Lines 57: computes `L` using `L=H+1i*R+1i*K`.
- Lines 60: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 63: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 66: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 67: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 68: computes `Cy` using `Cy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 71: computes `rho` using `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 74-75: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 79: computes `Dx` using `Dx=operator(spin_system,'Lx',parameters.decouple_f1{n})`.
- Lines 88: computes `rho_stack_cos_p` using `rho_stack_cos_p=step(spin_system,Hx,rho_stack,+pi/2)`.
- Lines 89: computes `rho_stack_sin_p` using `rho_stack_sin_p=step(spin_system,Hy,rho_stack,+pi/2)`.
- Lines 90: computes `rho_stack_cos_m` using `rho_stack_cos_m=step(spin_system,Hx,rho_stack,-pi/2)`.
- Lines 91: computes `rho_stack_sin_m` using `rho_stack_sin_m=step(spin_system,Hy,rho_stack,-pi/2)`.
- Lines 118: computes `rho_stack_cos` using `rho_stack_cos=rho_stack_cos_p-rho_stack_cos_m`.
- Lines 119: computes `rho_stack_sin` using `rho_stack_sin=rho_stack_sin_p-rho_stack_sin_m`.
- Lines 122: computes `[L,rho_stack_cos]` using `[L,rho_stack_cos]=decouple(spin_system,L,rho_stack_cos,parameters.spins(1))`.
- Lines 123: computes `[L,rho_stack_sin]` using `[L,rho_stack_sin]=decouple(spin_system,L,rho_stack_sin,parameters.spins(1))`.
- Lines 126: computes `fid.cos` using `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 132: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=hoesy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep two sweep widths, Hz
- parameters.npoints number of FID points for both
- dimensions
- parameters.spins nuclei on which the sequence runs,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei to decouple in F1, e.g.
- {'1H','13C'}
- parameters.tmix mixing time, seconds
- parameters.rho0 initial state
- parameters.needs should be set to {'rho_eq'}, this
- sequence needs the thermal equili-
- brium state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos,fid.sin -two components of the FID for F1 hyper-
- complex processing
- Note: this is an ideal heteronuclear NOESY model. Gradient and
- diffusion attenuation, finite-pulse losses, and experimental
- normalisation are outside this pulse sequence function.

## Implementation structure

- Phase-sensitive heteronuclear NOESY pulse sequence from:
- fid=hoesy(spin_system,parameters,H,R,K)
- parameters.sweep two sweep widths, Hz
- parameters.npoints number of FID points for both
- dimensions
- parameters.spins nuclei on which the sequence runs,
- e.g. {'15N','13C'}
- parameters.decouple_f1 nuclei to decouple in F1, e.g.
- {'1H','13C'}
- parameters.tmix mixing time, seconds
- parameters.rho0 initial state
- parameters.needs should be set to {'rho_eq'}, this

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `timestep()`, `homospoil()`, `decouple()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
