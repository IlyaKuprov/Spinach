# experiments/spen/psycosy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/spen/psycosy.m`
- Signature: `fid=psycosy(spin_system,parameters,H,R,K,G,F)`
- Total lines: 253

## Purpose

Alan Kenwright's spatially encoded COSY sequence described in fid=psycosy(spin_system,parameters,H,R,K,G,F)

## Physical / mathematical content

- SPEN experiment implementations. These files combine shaped pulses, gradients, spatial encoding, and often diffusion-aware propagation.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix mixing time, seconds
- parameters.gamp gradient amplitude, T/m
- parameters.sal_ang flip angle of the saltire chirp (degrees)
- parameters.sal_dur pulse width of saltire chirp (s)
- parameters.sal_del chirp pulse gradient duration (s)
- parameters.sal_swp sweep width of saltire chirp (Hz)
- parameters.sal_npt number of points in the saltire chirp
- parameters.sal_smf saltire chirp smoothing factor
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

- fid -two-dimensional free induction decay

## Implementation structure

- Alan Kenwright's spatially encoded COSY sequence described in
- fid=psycosy(spin_system,parameters,H,R,K,G,F)
- parameters.sweep sweep width in Hz
- parameters.npoints number of points for both dimensions
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.tmix mixing time, seconds
- parameters.gamp gradient amplitude, T/m
- parameters.sal_ang flip angle of the saltire chirp (degrees)
- parameters.sal_dur pulse width of saltire chirp (s)
- parameters.sal_del chirp pulse gradient duration (s)
- parameters.sal_swp sweep width of saltire chirp (Hz)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `operator()`, `speye()`, `chirp_pulse()`, `cosd()`, `step()`, `evolution()`, `coherence()`, `shaped_pulse_xy()`, `ismember()`, `ismatrix()`, `all()`, `iscell()`, `isfield()`, `elseif()`, `ischar()`.
