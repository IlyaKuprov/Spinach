# experiments/nmr_liquids/ecosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/ecosy.m`
- Signature: `fid=ecosy(spin_system,parameters,H,R,K)`
- Total lines: 128

## Purpose

Phase-sensitive E.COSY pulse sequence from:

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

- Lines 40-41: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 43-44: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 46-47: Compute evolution timestep; implemented by `timestep=1/parameters.sweep`.
- Lines 49-50: Initial (post-pulse) and detection states; implemented by `rho0=state(spin_system,'Lx',parameters.spins{1})`.
- Lines 53-54: Get pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 57-59: Run F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho0,timestep, parameters.npoints(1)-1,'trajectory')`.
- Lines 61-62: Apply the second pulse (States quadrature); implemented by `rho_stack_cos=step(spin_system,Lx,rho_stack,pi/2)`.
- Lines 65-70: Apply the multiple-quantum filter; implemented by `rho_stack_cos=1*coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+2,-2]}})+ 2*coherence(spin_system,rho_stack_cos,{{parameters.spins{1},[+3,-3]}})+ 4*coherence…`.
- Lines 77-78: Apply the third pulse; implemented by `rho_stack_cos=step(spin_system,Lx,rho_stack_cos,pi/2)`.
- Lines 81-83: Run the F2 evolution; implemented by `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep, parameters.npoints(2)-1,'observable')`.

### Key state/data transformations

- Lines 44: computes `L` using `L=H+1i*R+1i*K`.
- Lines 47: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 50: computes `rho0` using `rho0=state(spin_system,'Lx',parameters.spins{1})`.
- Lines 51: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1})`.
- Lines 54: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 55: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 58-59: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho0,timestep, parameters.npoints(1)-1,'trajectory')`.
- Lines 62: computes `rho_stack_cos` using `rho_stack_cos=step(spin_system,Lx,rho_stack,pi/2)`.
- Lines 63: computes `rho_stack_sin` using `rho_stack_sin=step(spin_system,Ly,rho_stack,pi/2)`.
- Lines 82-83: computes `fid.cos` using `fid.cos=evolution(spin_system,L,coil,rho_stack_cos,timestep, parameters.npoints(2)-1,'observable')`.
- Lines 84-85: computes `fid.sin` using `fid.sin=evolution(spin_system,L,coil,rho_stack_sin,timestep, parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 90: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=ecosy(spin_system,parameters,H,R,K)
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

- fid.cos, fid.sin -real and imaginary components of the
- States quadrature signal
- Notes: implemented up to six-quantum orders as per the original
- paper. Let us know if you need more.

## Implementation structure

- Phase-sensitive E.COSY pulse sequence from:
- fid=ecosy(spin_system,parameters,H,R,K)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid.cos, fid.sin -real and imaginary components of the
- States quadrature signal
- paper. Let us know if you need more.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `state()`, `operator()`, `evolution()`, `step()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `any()`.
