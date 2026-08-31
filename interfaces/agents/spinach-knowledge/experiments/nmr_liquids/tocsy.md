# experiments/nmr_liquids/tocsy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/tocsy.m`
- Signature: `fid=tocsy(spin_system,parameters,H,R,K)`
- Total lines: 147

## Purpose

Amplitude-mode homonuclear TOCSY pulse sequence from:

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

- Lines 48-49: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 51-52: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 54-55: Coherent evolution timestep; implemented by `timestep=1./parameters.sweep`.
- Lines 57-58: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 60-61: Pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 64-65: First pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 67-68: F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 70-71: Mixing time under spin-lock; implemented by `rho_stack_cos=evolution(spin_system,L+2*pi*parameters.lamp*Lx,[],rho_stack,parameters.tmix,1,'final')`.
- Lines 74-75: F2 evolution; implemented by `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 52: computes `L` using `L=H+1i*R+1i*K`.
- Lines 55: computes `timestep` using `timestep=1./parameters.sweep`.
- Lines 58: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 61: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 62: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 65: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 68: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 71: computes `rho_stack_cos` using `rho_stack_cos=evolution(spin_system,L+2*pi*parameters.lamp*Lx,[],rho_stack,parameters.tmix,1,'final')`.
- Lines 72: computes `rho_stack_sin` using `rho_stack_sin=evolution(spin_system,L+2*pi*parameters.lamp*Ly,[],rho_stack,parameters.tmix,1,'final')`.
- Lines 75: computes `fid.cos` using `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep(2),parameters.npoints(2)-1,'observable')`.
- Lines 76: computes `fid.sin` using `fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Syntax

```matlab
fid=tocsy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep -sweep widths, Hz
- parameters.npoints -number of points for both
- dimensions
- parameters.spins -nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix -mixing time, seconds
- parameters.lamp -spin-lock power, Hz
- parameters.rho0 -initial state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid.cos, fid.sin -sine and cosine components
- of the States quadrature
- Note: this is an ideal continuous-spin-lock TOCSY, not an
- explicit MLEV, DIPSI, WALTZ, or clean-TOCSY composite
- pulse-train simulation.

## Implementation structure

- Amplitude-mode homonuclear TOCSY pulse sequence from:
- fid=tocsy(spin_system,parameters,H,R,K)
- parameters.sweep - sweep widths, Hz
- parameters.npoints - number of points for both
- dimensions
- parameters.spins - nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix - mixing time, seconds
- parameters.lamp - spin-lock power, Hz
- parameters.rho0 - initial state
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `timestep()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
