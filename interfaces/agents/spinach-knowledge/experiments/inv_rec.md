# experiments/inv_rec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/inv_rec.m`
- Signature: `fids=inv_rec(spin_system,parameters,H,R,K)`
- Total lines: 128

## Purpose

Inversion-recovery pulse sequence. Syntax: fids=inv_rec(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.sweep spectrum sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.max_delay longest relaxation delay
- parameters.n_delays number of relaxation delays to run
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fids -free induction decays for each delay starting from zero,
- a matrix with individual FIDs in columns
- Note: the relaxation superoperator must be thermalised.
- Zak El-Machachi

## Implementation structure

- Inversion-recovery pulse sequence. Syntax:
- fids=inv_rec(spin_system,parameters,H,R,K)
- parameters.sweep spectrum sweep width, Hz
- parameters.npoints number of points in the FID
- parameters.spins nuclei on which the sequence runs,
- specified as {'1H'}, {'13C'}, etc.
- parameters.max_delay longest relaxation delay
- parameters.n_delays number of relaxation delays to run
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- fids -free induction decays for each delay starting from zero,

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `equilibrium()`, `state()`, `operator()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `iscell()`, `ischar()`, `ismember()`, `isscalar()`, `strcmp()`.
