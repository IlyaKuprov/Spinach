# experiments/nmr_liquids/inept.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/inept.m`
- Signature: `fid=inept(spin_system,parameters,H,R,K)`
- Total lines: 146

## Purpose

Non-refocused INEPT pulse sequence. This returns the directly acquired coupled antiphase spectrum rather than a refocused, broadband-decoupled INEPT variant. Implemented as here:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 44-45: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 47-48: Get timing parameters; implemented by `tau=abs(1/(4*parameters.J))`.
- Lines 51-52: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 54-55: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 57-58: Pulse operators; implemented by `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 63-64: 90x pulse on H; implemented by `rho=step(spin_system,Hx,rho,pi/2)`.
- Lines 66-67: tau evolution; implemented by `rho=evolution(spin_system,L,[],rho,tau,1,'final')`.
- Lines 69-70: Two inversion pulses; implemented by `rho=step(spin_system,Cy+Hy,rho,pi)`.
- Lines 72-73: Second tau evolution; implemented by `rho=evolution(spin_system,L,[],rho,tau,1,'final')`.
- Lines 75-76: 90x pulse on C; implemented by `rho=step(spin_system,Cx,rho,pi/2)`.
- Lines 78-79: Split phase 90y pulses on H; implemented by `rho_pos=step(spin_system,Hy,rho,+pi/2)`.
- Lines 82-83: Phase cycle; implemented by `rho=(rho_pos-rho_neg)/2`.
- Lines 85-87: Detection; implemented by `fid=evolution(spin_system,L,coil,rho,timestep, parameters.npoints-1,'observable')`.

### Key state/data transformations

- Lines 45: computes `L` using `L=H+1i*R+1i*K`.
- Lines 48: computes `tau` using `tau=abs(1/(4*parameters.J))`.
- Lines 49: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 52: computes `rho` using `rho=equilibrium(spin_system)`.
- Lines 55: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 58: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 59: computes `Cy` using `Cy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 60: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 61: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 79: computes `rho_pos` using `rho_pos=step(spin_system,Hy,rho,+pi/2)`.
- Lines 80: computes `rho_neg` using `rho_neg=step(spin_system,Hy,rho,-pi/2)`.
- Lines 86-87: computes `fid` using `fid=evolution(spin_system,L,coil,rho,timestep, parameters.npoints-1,'observable')`.

### Local helper functions

- Line 92: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Syntax

```matlab
fid=inept(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1 F2} working nuclei,
- e.g. {'15N','1H'}
- parameters.J working scalar coupling
- in Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay
- Note: use dilute() to generate carbon isotopomers.
- Andrew Porter, Ilya Kuprov

## Implementation structure

- Non-refocused INEPT pulse sequence. This returns the directly
- acquired coupled antiphase spectrum rather than a refocused,
- broadband-decoupled INEPT variant. Implemented as here:
- fid=inept(spin_system,parameters,H,R,K)
- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1 F2} working nuclei,
- e.g. {'15N','1H'}
- parameters.J working scalar coupling
- in Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `equilibrium()`, `state()`, `operator()`, `step()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `strcmp()`, `any()`.
