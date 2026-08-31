# experiments/nmr_liquids/hetcor.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/hetcor.m`
- Signature: `fid=hetcor(spin_system,parameters,H,R,K)`
- Total lines: 176

## Purpose

Magnitude-mode HETCOR pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 49-50: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 52-53: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 55-56: Evolution timesteps; implemented by `timestep=1./parameters.sweep`.
- Lines 58-59: J-coupling evolution times; implemented by `delta_2=abs(1/(2*parameters.J))`.
- Lines 62-63: Initial state; implemented by `rho=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 65-66: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{2})`.
- Lines 68-69: Pulse operators; implemented by `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 73-75: First pulse on F1 with coherence selection; implemented by `rho= step(spin_system,Lx_F1,rho,pi/2)+ 1i*step(spin_system,Ly_F1,rho,pi/2)`.
- Lines 77-79: F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 84-85: 1/2J part of delta evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,delta_2,1,'final')`.
- Lines 87-89: Two pi/2 pulses with coherence selection; implemented by `rho_stack=+step(spin_system,Lx_F1+Lx_F2,rho_stack,+pi/2) +step(spin_system,Lx_F1+Lx_F2,rho_stack,-pi/2)`.
- Lines 91-92: 1/3J part of delta evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,delta_3,1,'final')`.
- Lines 94-95: F1 decoupling; implemented by `[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.decouple)`.
- Lines 97-99: F2 detection; implemented by `fid=evolution(spin_system,L,coil,rho_stack,timestep(2), parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 53: computes `L` using `L=H+1i*R+1i*K`.
- Lines 56: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 59: computes `delta_2` using `delta_2=abs(1/(2*parameters.J))`.
- Lines 60: computes `delta_3` using `delta_3=abs(1/(3*parameters.J))`.
- Lines 63: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1})`.
- Lines 66: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{2})`.
- Lines 69: computes `Lx_F1` using `Lx_F1=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 70: computes `Ly_F1` using `Ly_F1=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 71: computes `Lx_F2` using `Lx_F2=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 78-79: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 95: computes `[L,rho_stack]` using `[L,rho_stack]=decouple(spin_system,L,rho_stack,parameters.decouple)`.
- Lines 98-99: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,timestep(2), parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 104: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=hetcor(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths in the
- two frequency directions, Hz
- parameters.npoints [F1 F2] numbers of points in
- the two time directions in fid
- parameters.spins {F1 F2} nuclei (e.g. {'1H','13C'})
- parameters.decouple list of nuclei that detection
- time decoupling should be applied
- to -cell array of strings, e.g.
- {'1H','15N'})
- parameters.J working scalar coupling, Hz
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay for magnitude-
- mode processing
- Note: natural abundance experiments should make use of the iso-
- tope dilution functionality. See dilute.m function.
- Note: the fixed transfer delays are delta_2=1/(2J) and
- delta_3=1/(3J), using the absolute value of J.

## Implementation structure

- Magnitude-mode HETCOR pulse sequence from:
- fid=hetcor(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths in the
- two frequency directions, Hz
- parameters.npoints [F1 F2] numbers of points in
- the two time directions in fid
- parameters.spins {F1 F2} nuclei (e.g. {'1H','13C'})
- parameters.decouple list of nuclei that detection
- time decoupling should be applied
- to -cell array of strings, e.g.
- {'1H','15N'})
- parameters.J working scalar coupling, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `timestep()`, `decouple()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`, `strcmp()`.
