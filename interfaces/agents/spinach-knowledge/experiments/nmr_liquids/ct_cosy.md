# experiments/nmr_liquids/ct_cosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/ct_cosy.m`
- Signature: `fid=ct_cosy(spin_system,parameters,H,R,K)`
- Total lines: 164

## Purpose

Constant-time COSY pulse sequence with analytical coherence selec- tion, as described in: The F1 trace is flipped at the end to preserve the conventional indirect-dimension sign for the reversed-delay implementation.

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Set default final pulse angle; implemented by `if ~isfield(parameters,'angle')`.
- Lines 50-51: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 53-54: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 56-57: Initial condition; implemented by `if ~isfield(parameters,'rho0')`.
- Lines 61-62: Detection state; implemented by `if ~isfield(parameters,'coil')`.
- Lines 66-67: Get the time grid for the CT period; implemented by `t1_grid=(0:(parameters.npoints(1)-1))/parameters.sweep(1)`.
- Lines 70-71: Get the pulse operator; implemented by `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 73-74: Apply the first 90 degree pulse; implemented by `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 76-77: Select "+1" coherence; implemented by `rho=coherence(spin_system,rho,{{parameters.spins{1},+1}})`.
- Lines 79-80: Preallocate rho stack; implemented by `rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i)`.
- Lines 82-83: Loop over the value of t1; implemented by `parfor n=1:parameters.npoints(1)`.
- Lines 85-86: Run the first delay; implemented by `rho_current=evolution(spin_system,L,[],rho,CT/2-t1_grid(n)/2,1,'final')`.
- Lines 88-89: Run the pi pulse; implemented by `rho_current=step(spin_system,Hx,rho_current,pi)`.
- Lines 91-92: Run the second delay; implemented by `rho_current=evolution(spin_system,L,[],rho_current,CT/2+t1_grid(n)/2,1,'final')`.
- Lines 94-95: Assign the stack element; implemented by `rho_stack(:,n)=rho_current`.
- Lines 99-100: Final pulse; implemented by `rho_stack=step(spin_system,Hx,rho_stack,parameters.angle)`.
- Lines 102-104: Run the F2 evolution; implemented by `fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.sweep(2), parameters.npoints(2)-1,'observable')`.
- Lines 106-107: Flip t1 direction; implemented by `fid=fliplr(fid)`.

### Control flow inferred from the code

- Line 46: conditional branch on `~isfield(parameters,'angle')`.
- Line 57: conditional branch on `~isfield(parameters,'rho0')`.
- Line 62: conditional branch on `~isfield(parameters,'coil')`.
- Line 83: `parfor` loop over `n=1:parameters.npoints(1)`.

### Key state/data transformations

- Lines 47: computes `parameters.angle` using `parameters.angle=pi/2`.
- Lines 54: computes `L` using `L=H+1i*R+1i*K`.
- Lines 58: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz',parameters.spins{1},'cheap')`.
- Lines 63: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 67: computes `t1_grid` using `t1_grid=(0:(parameters.npoints(1)-1))/parameters.sweep(1)`.
- Lines 68: computes `CT` using `CT=t1_grid(end)`.
- Lines 71: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 74: computes `rho` using `rho=step(spin_system,Hx,parameters.rho0,pi/2)`.
- Lines 80: computes `rho_stack` using `rho_stack=zeros([size(rho,1) parameters.npoints(1)],'like',1i)`.
- Lines 86: computes `rho_current` using `rho_current=evolution(spin_system,L,[],rho,CT/2-t1_grid(n)/2,1,'final')`.
- Lines 95: computes `rho_stack(:,n)` using `rho_stack(:,n)=rho_current`.
- Lines 103-104: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho_stack,1/parameters.sweep(2), parameters.npoints(2)-1,'observable')`.

### Local helper functions

- Line 112: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalism.')`.

## Syntax

```matlab
fid=ct_cosy(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep a two-element vector giving the
- sweep widths in F1 and F2
- parameters.npoints a two-element vector giving the
- number of points in F1 and F2
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle final pulse angle in radians, defaults
- to pi/2
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -two-dimensional free induction decay

## Implementation structure

- Constant-time COSY pulse sequence with analytical coherence selec-
- tion, as described in:
- The F1 trace is flipped at the end to preserve the conventional
- indirect-dimension sign for the reversed-delay implementation.
- fid=ct_cosy(spin_system,parameters,H,R,K)
- parameters.sweep a two-element vector giving the
- sweep widths in F1 and F2
- parameters.npoints a two-element vector giving the
- number of points in F1 and F2
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.angle final pulse angle in radians, defaults

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isfield()`, `grumble()`, `state()`, `t1_grid()`, `operator()`, `step()`, `coherence()`, `evolution()`, `rho_stack()`, `fliplr()`, `ismember()`, `ismatrix()`, `all()`, `elseif()`, `any()`, `iscell()`.
