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
