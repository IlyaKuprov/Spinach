# experiments/nmr_liquids/dept.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/dept.m`
- Signature: `fid=dept(spin_system,parameters,H,R,K)`
- Total lines: 184

## Purpose

DEPT pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 46-47: Consistency check; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 49-50: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 52-53: Get timing parameters; implemented by `tau_j=abs(1/(2*parameters.J))`.
- Lines 56-57: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 59-60: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 62-63: Pulse operators; implemented by `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 67-68: Proton pulse; implemented by `rho=step(spin_system,Hx,rho,pi/2)`.
- Lines 70-71: J-coupling evolution; implemented by `rho=step(spin_system,L,rho,tau_j)`.
- Lines 73-74: Proton and carbon pulses; implemented by `rho=step(spin_system,Cx,rho,pi/2)`.
- Lines 78-79: Second tau evolution; implemented by `rho_a=step(spin_system,L,rho_a,tau_j)`.
- Lines 82-84: Proton and carbon pulses; implemented by `rho=step(spin_system,Cx,rho_a,+pi)+ step(spin_system,Cx,rho_b,-pi)`.
- Lines 87-88: Third tau evolution; implemented by `rho=step(spin_system,L,rho,tau_j)`.
- Lines 90-91: Decoupling; implemented by `[L,rho]=decouple(spin_system,L,rho,parameters.spins(2))`.
- Lines 93-94: Detection; implemented by `fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable')`.

### Key state/data transformations

- Lines 50: computes `L` using `L=H+1i*R+1i*K`.
- Lines 53: computes `tau_j` using `tau_j=abs(1/(2*parameters.J))`.
- Lines 54: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 57: computes `rho` using `rho=equilibrium(spin_system)`.
- Lines 60: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 63: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 64: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 65: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 75: computes `rho_a` using `rho_a=step(spin_system,Hx,rho,+pi)-step(spin_system,Hy,rho,+pi)`.
- Lines 76: computes `rho_b` using `rho_b=step(spin_system,Hx,rho,-pi)-step(spin_system,Hy,rho,-pi)`.
- Lines 91: computes `[L,rho]` using `[L,rho]=decouple(spin_system,L,rho,parameters.spins(2))`.
- Lines 94: computes `fid` using `fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 99: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Syntax

```matlab
fid=dept(spin_system,parameters,H,R,K)
```

## Parameters / inputs

- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1,F2} nuclei, e.g. {'13C','1H'}
- parameters.J working J-coupling in Hz
- parameters.beta the angle used in the selection
- pulse, radians
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -free induction decay
- Note: DEPT135 yields spectra with CH and CH3 signals in opposite phase
- to CH2 signals; DEPT90 yields spectra with only CH signals; DEPT45
- yields spectra with positive CH, CH2, and CH3 signals; quaternary
- carbons do not appear.
- Note: use dilute.m to generate carbon isotopomers.

## Implementation structure

- DEPT pulse sequence from:
- fid=dept(spin_system,parameters,H,R,K)
- parameters.sweep [F1] Sweep width in Hz
- parameters.npoints [F1] number of points
- parameters.spins {F1,F2} nuclei, e.g. {'13C','1H'}
- parameters.J working J-coupling in Hz
- parameters.beta the angle used in the selection
- pulse, radians
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fid -free induction decay

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `equilibrium()`, `state()`, `operator()`, `step()`, `decouple()`, `evolution()`, `ismember()`, `ismatrix()`, `all()`, `isfield()`, `elseif()`, `iscell()`, `ischar()`, `strcmp()`, `any()`.
