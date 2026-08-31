# experiments/bruker/hsqcetgpsi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/bruker/hsqcetgpsi.m`
- Signature: `fid=hsqcetgpsi(spin_system,parameters,H,R,K)`
- Total lines: 294

## Purpose

Sensitivity-improved echo/antiecho gradient-selected HSQC pulse sequence, based on the Bruker hsqcetgpsi pulse program and the standard HSQC sequence from: The gradient selection is represented analytically by coherence order selection statements

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 59-60: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 62-63: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 65-66: Coherent evolution timesteps; implemented by `timestep=1./parameters.sweep`.
- Lines 68-69: J-coupling evolution time; implemented by `delta=abs(1/(2*parameters.J))`.
- Lines 71-72: Initial condition; implemented by `rho=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 74-75: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 77-78: Pulse operators; implemented by `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 83-84: Initial proton pulse; implemented by `rho=step(spin_system,Hx,rho,pi/2)`.
- Lines 86-87: First INEPT evolution; implemented by `rho=evolution(spin_system,L,[],rho,delta/2,1,'final')`.
- Lines 89-90: Refocusing pulses; implemented by `rho=step(spin_system,Hx+Cx,rho,pi)`.
- Lines 92-93: Second INEPT evolution; implemented by `rho=evolution(spin_system,L,[],rho,delta/2,1,'final')`.
- Lines 95-96: Proton trim pulse; implemented by `rho=step(spin_system,Hx,rho,parameters.trim_angle)`.
- Lines 98-99: Transfer pulses; implemented by `rho=step(spin_system,Hy,rho,pi/2)`.
- Lines 102-104: First half of F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 106-107: F1 dimension refocusing pulses; implemented by `for n=1:numel(parameters.decouple_f1)`.
- Lines 112-114: Second half of F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,timestep(1)/2, parameters.npoints(1)-1,'refocus')`.
- Lines 116-118: First analytical gradient selection; implemented by `rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2},0}, {parameters.spins{1},+1}})`.
- Lines 122-123: Carbon inversion pulse after the first gradient; implemented by `rho_stack_pos=step(spin_system,Cx,rho_stack_pos,pi)`.

### Control flow inferred from the code

- Line 107: `for` loop over `n=1:numel(parameters.decouple_f1)`.

### Key state/data transformations

- Lines 63: computes `L` using `L=H+1i*R+1i*K`.
- Lines 66: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 69: computes `delta` using `delta=abs(1/(2*parameters.J))`.
- Lines 72: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{2},'cheap')`.
- Lines 75: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{2},'cheap')`.
- Lines 78: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 79: computes `Cy` using `Cy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 80: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 81: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 103-104: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1)/2, parameters.npoints(1)-1,'trajectory')`.
- Lines 108: computes `Lp` using `Lp=operator(spin_system,'L+',parameters.decouple_f1{n})`.
- Lines 117-118: computes `rho_stack_pos` using `rho_stack_pos=coherence(spin_system,rho_stack,{{parameters.spins{2},0}, {parameters.spins{1},+1}})`.
- Lines 119-120: computes `rho_stack_neg` using `rho_stack_neg=coherence(spin_system,rho_stack,{{parameters.spins{2},0}, {parameters.spins{1},-1}})`.
- Lines 181: computes `[L,rho_stack_pos]` using `[L,rho_stack_pos]=decouple(spin_system,L,rho_stack_pos,parameters.decouple_f2)`.
- Lines 182: computes `[L,rho_stack_neg]` using `[L,rho_stack_neg]=decouple(spin_system,L,rho_stack_neg,parameters.decouple_f2)`.
- Lines 185-186: computes `fid.pos` using `fid.pos=evolution(spin_system,L,coil,rho_stack_pos, timestep(2),parameters.npoints(2)-1,'observable')`.
- Lines 187-188: computes `fid.neg` using `fid.neg=evolution(spin_system,L,coil,rho_stack_neg, timestep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 193: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=hsqcetgpsi(spin_system,parameters,H,R,K)
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
- parameters.trim_angle proton trim pulse angle, rad
- parameters.si_time sensitivity improvement delay, s
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.pos,fid.neg -echo and antiecho components of the
- signal.
- Note: natural abundance simulations should make use of the isotope
- dilution functionality. See dilute.m function.

## Implementation structure

- Sensitivity-improved echo/antiecho gradient-selected HSQC pulse
- sequence, based on the Bruker hsqcetgpsi pulse program and the
- standard HSQC sequence from:
- The gradient selection is represented analytically by coherence order
- selection statements
- fid=hsqcetgpsi(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths, Hz
- parameters.npoints [F1 F2] numbers of points
- parameters.spins {F1 F2} nuclei (e.g. '13C','1H')
- parameters.decouple_f2 nuclei to decouple in F2, e.g.
- {'15N','13C'}
- parameters.decouple_f1 nuclei that receive midpoint

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `timestep()`, `coherence()`, `decouple()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
