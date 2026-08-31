# experiments/nmr_liquids/inadequate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/inadequate.m`
- Signature: `fid=inadequate(spin_system,parameters,H,R,K)`
- Total lines: 158

## Purpose

INADEQUATE pulse sequence. Selects double-quantum coherence from coupled carbon pairs, then converts it back into observable single- quantum magnetisation. At natural abundance 13C, this produces only 13C pair subspectra. Implemented as described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 46-47: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 49-50: Decoupling; implemented by `L=decouple(spin_system,L,[],parameters.decouple)`.
- Lines 52-53: Sequence timing; implemented by `timestep=1./parameters.sweep`.
- Lines 56-57: Initial and detection states; implemented by `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 60-61: Pulse operators; implemented by `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 64-65: Pulse 90x; implemented by `rho=step(spin_system,Cx,rho,pi/2)`.
- Lines 67-68: J-coupling evolution; implemented by `rho=step(spin_system,L,rho,tau)`.
- Lines 70-71: Pulse 180y; implemented by `rho=step(spin_system,Cy,rho,pi)`.
- Lines 79-80: Select double-quantum coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},[2 -2]}})`.
- Lines 82-83: Pulse on 90x; implemented by `rho=step(spin_system,Cx,rho,pi/2)`.
- Lines 85-87: Detection; implemented by `fid=evolution(spin_system,L,coil,rho,timestep(1), parameters.npoints(1)-1,'observable')`.

### Key state/data transformations

- Lines 47: computes `L` using `L=H+1i*R+1i*K`.
- Lines 53: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 54: computes `tau` using `tau=abs(1/(4*parameters.J))`.
- Lines 57: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 58: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 61: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 62: computes `Cy` using `Cy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 86-87: computes `fid` using `fid=evolution(spin_system,L,coil,rho,timestep(1), parameters.npoints(1)-1,'observable')`.

### Local helper functions

- Line 92: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Syntax

```matlab
fid=inadequate(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points in the fid
- parameters.spins active nuclei, e.g. {'13C'}
- parameters.decouple nuclei to decouple, e.g. {'1H'}
- parameters.J working J-coupling in Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- fid -free induction decay
- Note: use dilute.m to generate carbon pair isotopomers.

## Implementation structure

- INADEQUATE pulse sequence. Selects double-quantum coherence from
- coupled carbon pairs, then converts it back into observable single-
- quantum magnetisation. At natural abundance 13C, this produces only
- 13C pair subspectra. Implemented as described in:
- fid=inadequate(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points in the fid
- parameters.spins active nuclei, e.g. {'13C'}
- parameters.decouple nuclei to decouple, e.g. {'1H'}
- parameters.J working J-coupling in Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `state()`, `operator()`, `step()`, `coherence()`, `evolution()`, `timestep()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`.
