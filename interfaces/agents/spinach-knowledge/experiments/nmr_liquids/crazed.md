# experiments/nmr_liquids/crazed.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/crazed.m`
- Signature: `fid=crazed(spin_system,parameters,H,R,K)`
- Total lines: 138

## Purpose

CRAZED pulse sequence. Ideal analytical coherence-pathway version of the sequence described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 47-48: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 50-51: Compute the evolution timestep; implemented by `timestep=1/parameters.sweep`.
- Lines 53-54: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1})`.
- Lines 56-57: Get the pulse operator; implemented by `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 59-60: Apply the first pulse; implemented by `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 62-64: Run the F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,timestep, parameters.npoints(1)-1,'trajectory')`.
- Lines 66-67: Apply a double-quantum filter:; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},2}})`.
- Lines 69-70: Apply the second pulse; implemented by `rho_stack=step(spin_system,Ly,rho_stack,parameters.angle)`.
- Lines 72-73: Apply a single-quantum filter; implemented by `rho_stack=coherence(spin_system,rho_stack,{{parameters.spins{1},1}})`.
- Lines 75-77: Run the F2 evolution; implemented by `fid=evolution(spin_system,L,coil,rho_stack,timestep, parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 48: computes `L` using `L=H+1i*R+1i*K`.
- Lines 51: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 54: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1})`.
- Lines 57: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 60: computes `rho` using `rho=step(spin_system,Ly,parameters.rho0,pi/2)`.
- Lines 63-64: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,timestep, parameters.npoints(1)-1,'trajectory')`.
- Lines 76-77: computes `fid` using `fid=evolution(spin_system,L,coil,rho_stack,timestep, parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 82: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=crazed(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle
- parameters.rho0 initial condition
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay
- Note: the gradient selection is represented by explicit coherence
- projection onto the +2 double-quantum branch followed by the
- +1 observable single-quantum branch.

## Implementation structure

- CRAZED pulse sequence. Ideal analytical coherence-pathway version
- of the sequence described in:
- fid=crazed(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle second pulse angle
- parameters.rho0 initial condition
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `any()`.
