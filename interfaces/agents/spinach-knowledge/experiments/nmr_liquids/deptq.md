# experiments/nmr_liquids/deptq.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/nmr_liquids/deptq.m`
- Signature: `fid=deptq(spin_system,parameters,H,R,K)`
- Total lines: 182

## Purpose

DEPTQ pulse sequence from:

## Physical / mathematical content

- Liquid-state pulse sequence implementations. These are production experiment kernels that carry out coherence transfer, mixing, refocusing, decoupling, and indirect evolution on precomputed Hamiltonian/relaxation/kinetics operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 48-49: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 51-52: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 54-55: Get timing parameters; implemented by `tau_j=abs(1/(2*parameters.J))`.
- Lines 58-59: Isotropic thermal equilibrium; implemented by `rho=equilibrium(spin_system)`.
- Lines 61-62: Detection state; implemented by `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 64-65: Pulse operators; implemented by `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 70-71: Carbon pulse; implemented by `rho=step(spin_system,Cy,rho,-pi/2)`.
- Lines 73-74: J-coupling evolution; implemented by `rho=step(spin_system,L,rho,tau_j)`.
- Lines 76-77: Proton and carbon pulses; implemented by `rho=step(spin_system,Hx,rho,+pi/2)`.
- Lines 83-84: J-coupling evolution; implemented by `rho_a=step(spin_system,L,rho_a,tau_j)`.
- Lines 89-91: Proton and carbon pulses; implemented by `rho_a=+step(spin_system,Hx,rho_a,+pi) -step(spin_system,Hy,rho_c,+pi)`.
- Lines 101-103: Proton and carbon pulses; implemented by `rho=step(spin_system,Cx,rho_a,+pi)+ step(spin_system,Cx,rho_b,-pi)`.
- Lines 109-110: Detection pulse; implemented by `rho=step(spin_system,Cx,rho,-pi/2)`.
- Lines 112-113: Proton decoupling; implemented by `[L,rho]=decouple(spin_system,L,rho,parameters.spins(2))`.
- Lines 115-116: Phase cycled detection; implemented by `fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable')`.

### Key state/data transformations

- Lines 52: computes `L` using `L=H+1i*R+1i*K`.
- Lines 55: computes `tau_j` using `tau_j=abs(1/(2*parameters.J))`.
- Lines 56: computes `timestep` using `timestep=1/parameters.sweep`.
- Lines 59: computes `rho` using `rho=equilibrium(spin_system)`.
- Lines 62: computes `coil` using `coil=state(spin_system,'L+',parameters.spins{1},'cheap')`.
- Lines 65: computes `Cx` using `Cx=operator(spin_system,'Lx',parameters.spins{1})`.
- Lines 66: computes `Cy` using `Cy=operator(spin_system,'Ly',parameters.spins{1})`.
- Lines 67: computes `Hx` using `Hx=operator(spin_system,'Lx',parameters.spins{2})`.
- Lines 68: computes `Hy` using `Hy=operator(spin_system,'Ly',parameters.spins{2})`.
- Lines 78: computes `rho_a` using `rho_a=step(spin_system,Cx,rho,+pi)`.
- Lines 79: computes `rho_b` using `rho_b=step(spin_system,Cx,rho,-pi)`.
- Lines 80: computes `rho_c` using `rho_c=step(spin_system,Cy,rho,+pi)`.
- Lines 81: computes `rho_d` using `rho_d=step(spin_system,Cy,rho,-pi)`.
- Lines 113: computes `[L,rho]` using `[L,rho]=decouple(spin_system,L,rho,parameters.spins(2))`.
- Lines 116: computes `fid` using `fid=evolution(spin_system,L,coil,rho,timestep,parameters.npoints-1,'observable')`.

### Local helper functions

- Line 121: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv formalisms.')`.

## Syntax

```matlab
fid=deptq(spin_system,parameters,H,R,K)
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
- Note: this implementation is the DEPTQ135-style fixed first
- proton pulse variant. The beta parameter controls the
- last proton editing pulse.
- Note: use dilute.m to generate carbon isotopomers.
- Note: the sequence differs from dept.m in that quaternary carbons
- do appear.

## Implementation structure

- DEPTQ pulse sequence from:
- fid=deptq(spin_system,parameters,H,R,K)
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
