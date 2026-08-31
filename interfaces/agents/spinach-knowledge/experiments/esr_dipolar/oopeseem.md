# experiments/esr_dipolar/oopeseem.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/oopeseem.m`
- Signature: `fid=oopeseem(spin_system,parameters,H,R,K)`
- Total lines: 125

## Purpose

Out-of-phase ESEEM pulse sequence with the first pulse set to pi/4 to probe two-electron correlations in the initial condition. Syntax: fid=eseem(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 45-46: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 48-49: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 51-52: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 54-55: Set the defaults; implemented by `if ~isfield(parameters,'screen'), parameters.screen=[]; end`.
- Lines 57-58: First pulse; implemented by `rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/4)`.
- Lines 60-62: Spin echo; implemented by `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep/2, (parameters.npoints-1),'trajectory',parameters.screen)`.
- Lines 63-64: Second pulse; implemented by `rho_stack=step(spin_system,parameters.pulse_op,rho_stack,pi)`.
- Lines 66-68: Spin echo; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep/2, (parameters.npoints-1),'refocus',parameters.coil)`.
- Lines 69-70: Detect; implemented by `fid=transpose(parameters.coil'*rho_stack)`.

### Control flow inferred from the code

- Line 55: conditional branch on `~isfield(parameters,'screen'), parameters.screen=[]; end`.

### Key state/data transformations

- Lines 46: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 52: computes `L` using `L=H+1i*R+1i*K`.
- Lines 58: computes `rho` using `rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/4)`.
- Lines 61-62: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep/2, (parameters.npoints-1),'trajectory',parameters.screen)`.
- Lines 70: computes `fid` using `fid=transpose(parameters.coil'*rho_stack)`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.npoints number of time points to be
- computed
- parameters.timestep simulation time step, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.screen optional screen state (must be
- the Hermitian conjugate of the
- detection state)
- parameters.pulse_op pulse operator A, the propagators
- for the pulses will be exp(-i*A*pi)
- and exp(-i*A*pi/4)
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -OOP-ESEEM time trace
- Note: the sequence uses ideal pulses, replace with shaped_pulse_af()
- to have soft pulses instead.

## Implementation structure

- Out-of-phase ESEEM pulse sequence with the first pulse set to pi/4 to
- probe two-electron correlations in the initial condition. Syntax:
- fid=eseem(spin_system,parameters,H,R,K)
- parameters.npoints number of time points to be
- computed
- parameters.timestep simulation time step, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.screen optional screen state (must be
- the Hermitian conjugate of the
- detection state)
- parameters.pulse_op pulse operator A, the propagators

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `isfield()`, `step()`, `evolution()`, `transpose()`, `ismember()`, `ismatrix()`, `all()`, `any()`.
