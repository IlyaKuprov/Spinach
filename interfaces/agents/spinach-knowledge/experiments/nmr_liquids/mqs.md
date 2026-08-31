# experiments/nmr_liquids/mqs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/mqs.m`
- Signature: `fid=mqs(spin_system,parameters,H,R,K)`
- Total lines: 178

## Purpose

2D multiple-quantum NMR pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 48-49: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 51-52: Get pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 55-56: Apply the first 90 deg pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 58-59: Apply first evolution period; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delay,1,'final')`.
- Lines 61-62: Apply 180 deg pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 64-65: Apply second evolution period; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delay,1,'final')`.
- Lines 67-68: Select operator; implemented by `if mod(parameters.mqorder,2)==0`.
- Lines 70-71: Apply the second 90 deg pulse about x; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 75-76: Apply the second 90 deg pulse about y; implemented by `rho=step(spin_system,Ly,rho,pi/2)`.
- Lines 80-81: Coherence selection; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},parameters.mqorder}})`.
- Lines 83-84: Run the F1 evolution; implemented by `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 86-87: Apply third 90 deg pulse; implemented by `rho_stack=step(spin_system,Lx,rho_stack,parameters.angle)`.
- Lines 89-90: Run the F2 evolution; implemented by `fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.sweep(2),parameters.npoints(2)-1,'observable')`.

### Control flow inferred from the code

- Line 68: conditional branch on `mod(parameters.mqorder,2)==0`.

### Key state/data transformations

- Lines 49: computes `L` using `L=H+1i*R+1i*K`.
- Lines 52: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 53: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 56: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 84: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,1/parameters.sweep(1),parameters.npoints(1)-1,'trajectory')`.
- Lines 90: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.sweep(2),parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 95: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=mqs(spin_system,parameters,H,R,K)
This function should be invoked through liquid.m context,
which would provide H, R, and K.
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of points
- parameters.spins working spins, e.g. {'1H','1H'}
- parameters.angle flip angle, radians
- parameters.delay J-coupling evolution delay, seconds
- parameters.mqorder coherence order to select
- parameters.rho0 initial state
- parameters.coil detection state

## Outputs

- fid -2D magnitude-mode free induction decay
- Note: this is the non-refocused multiple-quantum/MaxQ variant.
- Use mqs_refocus.m when post-mixing refocusing is required.

## Implementation structure

- 2D multiple-quantum NMR pulse sequence from:
- fid=mqs(spin_system,parameters,H,R,K)
- This function should be invoked through liquid.m context,
- which would provide H, R, and K.
- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of points
- parameters.spins working spins, e.g. {'1H','1H'}
- parameters.angle flip angle, radians
- parameters.delay J-coupling evolution delay, seconds
- parameters.mqorder coherence order to select
- parameters.rho0 initial state
- parameters.coil detection state

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`, `isscalar()`.
