# experiments/spen/psyche.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/psyche.m`
- Signature: `fid=psyche(spin_system,parameters,H,R,K,G,F)`
- Total lines: 266

## Purpose

PSYCHE pure-shift NMR pulse sequence. Syntax: fid=psyche_1d(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.rho0 initial state
- parameters.coil detection state
- parameters.spins nuclei on which the sequence runs
- parameters.g_amp gradient amplitude (T/m)
- parameters.dims size of the sample (m)
- parameters.npts number of discretization points in the grid
- parameters.sweep spectral range (Hz)
- parameters.npoints number of points in the sweep
- parameters.zerofill number of points for the zero filling
- parameters.diff diffusion constant (m^2/s)
- H Fokker-Planck Hamiltonian, received
- from the imaging context
- R Fokker-Planck relaxation superoperator,
- received from the imaging context
- K Fokker-Planck kinetics superoperator,
- received from the imaging context
- G Fokker-Planck gradient superoperators,
- received from the imaging context
- F Fokker-Planck diffusion and flow super-
- operator, received from the context

## Outputs

- fid -a PSYCHE free induction decay as a 2D array

## Implementation structure

- PSYCHE pure-shift NMR pulse sequence. Syntax:
- fid=psyche_1d(spin_system,parameters,H,R,K,G,F)
- parameters.rho0 initial state
- parameters.coil detection state
- parameters.spins nuclei on which the sequence runs
- parameters.g_amp gradient amplitude (T/m)
- parameters.dims size of the sample (m)
- parameters.npts number of discretization points in the grid
- parameters.sweep spectral range (Hz)
- parameters.npoints number of points in the sweep
- parameters.zerofill number of points for the zero filling
- parameters.diff diffusion constant (m^2/s)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `strcmp()`, `cosd()`, `step()`, `evolution()`, `coherence()`, `shaped_pulse_xy()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `elseif()`.
