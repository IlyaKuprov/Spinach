# experiments/nmr_liquids/mqs_refocus.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/mqs_refocus.m`
- Signature: `fid=mqs_refocus(spin_system,parameters,H,R,K)`
- Total lines: 205

## Purpose

Multiple quantum correlation pulse sequence with refocusing, as described in:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 54-55: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 57-58: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 60-61: Get pulse operators; implemented by `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 64-65: Apply the first 90 deg pulse; implemented by `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 67-68: Apply first evolution period; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delay_1,1,'final')`.
- Lines 70-71: Apply 180 deg pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 73-74: Apply second evolution period; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delay_1,1,'final')`.
- Lines 76-77: Select operator; implemented by `if mod(parameters.mqorder(1),2)==0`.
- Lines 79-80: Apply the second 90 deg pulse about x; implemented by `rho=step(spin_system,Lx,rho,pi/2)`.
- Lines 83-84: Apply the second 90 deg pulse about y; implemented by `rho=step(spin_system,Ly,rho,pi/2)`.
- Lines 87-88: Coherence selection; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},parameters.mqorder(1)}})`.
- Lines 90-92: Run the F1 evolution; implemented by `rho=evolution(spin_system,L,[],rho,1/parameters.sweep(1), parameters.npoints(1)-1,'trajectory')`.
- Lines 94-95: Apply third pulse; implemented by `rho=step(spin_system,Lx,rho,parameters.angle)`.
- Lines 97-98: Select coherence -1; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},parameters.mqorder(2)}})`.
- Lines 100-101: Apply first evolution period; implemented by `rho=evolution(spin_system,L,[],rho,parameters.delay_2,1,'final')`.
- Lines 103-104: Apply refocussing pulse; implemented by `rho=step(spin_system,Lx,rho,pi)`.
- Lines 109-111: Run the F2 evolution; implemented by `fid=evolution(spin_system,L,parameters.coil,rho,1/parameters.sweep(2), parameters.npoints(2)-1,'observable')`.

### Control flow inferred from the code

- Line 77: conditional branch on `mod(parameters.mqorder(1),2)==0`.

### Key state/data transformations

- Lines 58: computes `L` using `L=H+1i*R+1i*K`.
- Lines 61: computes `Lx` using `Lx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 62: computes `Ly` using `Ly=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 65: computes `rho` using `rho=step(spin_system,Lx,parameters.rho0,pi/2)`.
- Lines 110-111: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho,1/parameters.sweep(2), parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 116: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=mqs_refocus(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of fid points
- parameters.spins {F1 F2} nuclei, in this case:
- {'1H','1H'}
- parameters.angle flip angle for the final
- pulse, radians
- parameters.mqorder coherence orders to select,
- a two-element integer array
- parameters.delay_1 first evolution delay, seconds
- parameters.delay_2 second evolution delay, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- H -Hamiltonian superoperator, provided by the context function
- R -relaxation superoperator, provided by the context function
- K -kinetics superoperator, provided by the context function

## Outputs

- fid -2D free induction decay for amplitude-mode processing
- Note: this implementation is homonuclear in practice and uses
- exact analytical coherence-order projection rather than an
- explicit phase cycle.

## Implementation structure

- Multiple quantum correlation pulse sequence with refocusing,
- as described in:
- fid=mqs_refocus(spin_system,parameters,H,R,K)
- parameters.sweep [F1 F2] sweep widths (Hz)
- parameters.npoints [F1 F2] numbers of fid points
- parameters.spins {F1 F2} nuclei, in this case:
- {'1H','1H'}
- parameters.angle flip angle for the final
- pulse, radians
- parameters.mqorder coherence orders to select,
- a two-element integer array
- parameters.delay_1 first evolution delay, seconds

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `step()`, `evolution()`, `coherence()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `any()`, `iscell()`, `ischar()`.
