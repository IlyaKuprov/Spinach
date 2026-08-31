# experiments/esr_dipolar/eseem.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/eseem.m`
- Signature: `fid=eseem(spin_system,parameters,H,R,K)`
- Total lines: 120

## Purpose

ESEEM pulse sequence with ideal hard pulses. Syntax: fid=eseem(spin_system,parameters,H,R,K)

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 42-43: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 45-46: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 48-49: Set the defaults; implemented by `if ~isfield(parameters,'screen'), parameters.screen=[]; end`.
- Lines 51-52: First pulse; implemented by `rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/2)`.
- Lines 54-56: Spin echo; implemented by `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep/2, (parameters.npoints-1),'trajectory',parameters.screen)`.
- Lines 57-58: Second pulse; implemented by `rho_stack=step(spin_system,parameters.pulse_op,rho_stack,pi)`.
- Lines 60-62: Spin echo; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep/2, (parameters.npoints-1),'refocus',parameters.coil)`.
- Lines 63-64: Detect; implemented by `fid=transpose(full(parameters.coil'*rho_stack))`.

### Control flow inferred from the code

- Line 49: conditional branch on `~isfield(parameters,'screen'), parameters.screen=[]; end`.

### Key state/data transformations

- Lines 40: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 46: computes `L` using `L=H+1i*R+1i*K`.
- Lines 52: computes `rho` using `rho=step(spin_system,parameters.pulse_op,parameters.rho0,pi/2)`.
- Lines 55-56: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep/2, (parameters.npoints-1),'trajectory',parameters.screen)`.
- Lines 64: computes `fid` using `fid=transpose(full(parameters.coil'*rho_stack))`.

### Local helper functions

- Line 69: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})`.
  - Representative operation: `error('this function is only available for sphten-liouv and zeeman-liouv formalisms.')`.

## Parameters / inputs

- parameters.npoints number of points to be computed
- parameters.timestep simulation time step, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.screen optional screen state (must be
- the Hermitian conjugate of the
- detection state)
- parameters.pulse_op pulse operator
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fid -time domaing signal whose Fourier transform is the
- ESEEM spectrum

## Implementation structure

- ESEEM pulse sequence with ideal hard pulses. Syntax:
- fid=eseem(spin_system,parameters,H,R,K)
- parameters.npoints number of points to be computed
- parameters.timestep simulation time step, seconds
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.screen optional screen state (must be
- the Hermitian conjugate of the
- detection state)
- parameters.pulse_op pulse operator
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `isfield()`, `step()`, `evolution()`, `transpose()`, `ismember()`, `ismatrix()`, `all()`, `any()`.
