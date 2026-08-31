# experiments/esr_dipolar/sifter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/esr_dipolar/sifter.m`
- Signature: `fid=sifter(spin_system,parameters,H,R,K)`
- Total lines: 113

## Purpose

SIFTER pulse sequence. Syntax: fid=sifter(spin_system,parameters,H,R,K) where H is the Hamiltonian matrix, R is the relaxation matrix and K is the chemical kinetics matrix.

## Physical / mathematical content

- Dipolar ESR experiment implementations. The pulse logic resolves dipolar couplings by echo modulation, with selective excitation and time-domain accumulation.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(parameters,H,R,K)`.
- Lines 36-37: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 39-40: First 90-degree pulse along X; implemented by `rho=step(spin_system,parameters.pulse_opx,parameters.rho0,pi/2)`.
- Lines 42-44: First part of the t1 period; implemented by `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep, parameters.npoints/2-1,'trajectory')`.
- Lines 46-47: Second 180-degree pulse along +X; implemented by `rho_stack=step(spin_system,parameters.pulse_opx,rho_stack,pi)`.
- Lines 49-51: Second part of the t1 period; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack,parameters.timestep, parameters.npoints/2-1,'refocus')`.
- Lines 53-54: Third 90-degree pulse along +Y; implemented by `rho_stack=step(spin_system,parameters.pulse_opy,rho_stack,pi/2)`.
- Lines 56-58: First part of the t2 period; implemented by `rho_stack=evolution(spin_system,L,[],rho_stack(:,end:-1:1),parameters.timestep, parameters.npoints/2-1,'refocus')`.
- Lines 60-61: Fourth 180-degree pulse along +X; implemented by `rho_stack=step(spin_system,parameters.pulse_opx,rho_stack,pi)`.
- Lines 63-65: Second part of the t2 period, 2D detection mode; implemented by `fid=evolution(spin_system,L,parameters.coil,rho_stack,parameters.timestep, parameters.npoints/2-1,'observable')`.

### Key state/data transformations

- Lines 37: computes `L` using `L=H+1i*R+1i*K`.
- Lines 40: computes `rho` using `rho=step(spin_system,parameters.pulse_opx,parameters.rho0,pi/2)`.
- Lines 43-44: computes `rho_stack` using `rho_stack=evolution(spin_system,L,[],rho,parameters.timestep, parameters.npoints/2-1,'trajectory')`.
- Lines 64-65: computes `fid` using `fid=evolution(spin_system,L,parameters.coil,rho_stack,parameters.timestep, parameters.npoints/2-1,'observable')`.

### Local helper functions

- Line 70: `grumble()` — `function grumble(parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.npoints -number of points in time evolution
- parameters.timestep -simulation time step, seconds
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_opx -pulse operator in X phase
- parameters.pulse_opy -pulse operator in Y phase

## Outputs

- fid -a 2D free induction decay

## Implementation structure

- SIFTER pulse sequence. Syntax:
- fid=sifter(spin_system,parameters,H,R,K)
- where H is the Hamiltonian matrix, R is the relaxation matrix
- and K is the chemical kinetics matrix.
- parameters.npoints -number of points in time evolution
- parameters.timestep -simulation time step, seconds
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.pulse_opx -pulse operator in X phase
- parameters.pulse_opy -pulse operator in Y phase
- fid -a 2D free induction decay
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `step()`, `evolution()`, `rho_stack()`, `ismatrix()`, `all()`, `isfield()`, `isscalar()`.
