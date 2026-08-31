# experiments/nmr_liquids/cosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/cosy.m`
- Signature: `fid=cosy(spin_system,parameters,H,R,K)`
- Total lines: 132

## Purpose

Phase-sensitive COSY pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 47-48: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 50-51: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 53-54: Initial state up to an overall multiplier; implemented by `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 56-57: Detection state up to an overall multiplier; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 59-60: Get the pulse operator; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 62-63: Apply the first pulse; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 65-67: Run the F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep, parameters.npoints(1)-1,'trajectory')`.
- Lines 69-70: Select "+1" coherence; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},+1}})`.
- Lines 72-73: Apply the second pulse; implemented by `rho_stack=step(spin_system,Lx,rho_stack,parameters.angle)`.
- Lines 75-77: Run the F2 evolution; implemented by `fid=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep, parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 51: computes `L` using `L=H+1i*R+1i*K`.
- Lines 54: computes `rho` using `rho=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 57: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 60: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 66-67: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep, parameters.npoints(1)-1,'trajectory')`.
- Lines 76-77: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,1/parameters.sweep, parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 82: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle in radians, usu-
- ally pi/2, but also allows COSY45,
- COSY60, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay
- Note: if the second pulse angle differs from 90 degrees, magni-
- tude mode plotting is advised.
- Note: the implementation analytically retains the +1 t1 coher-
- ence order, corresponding to a single phase-sensitive
- pathway.

## Implementation structure

- Phase-sensitive COSY pulse sequence from:
- fid=cosy(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle in radians, usu-
- ally pi/2, but also allows COSY45,
- COSY60, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `any()`.
