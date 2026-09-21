# experiments/cp_acquire_soft.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/cp_acquire_soft.m`
- Signature: `fid=cp_acquire_soft(spin_system,parameters,H,R,K)`
- Total lines: 159

## Purpose

Cross-polarisation experiment in the rotating frame, followed by time-domain FID acquisition. The CP stage is preceded by wiping of the low-gamma spins and followed by FID acquisition with deco- upling of the high-gamma spins. Syntax: fid=cp_acquire_soft(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.spins -working spins, a cell array of
- strings with high-gamma spin fi-
- rst, and low-gamma spin last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- channel, Hz
- parameters.cp_pwr -nutation frequencies on the two
- channels during the CP contact
- time, a two-element vector, Hz
- parameters.cp_dur -duration of the contact time, s
- parameters.rho0 -initial state, the state of the
- low-gamma spins will be wiped
- parameters.coil -detection state
- parameters.sweep -sweep width for the FID, Hz
- parameters.npoints -number of points in the FID
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- fid -signal detected on the coil state during
- system evolution

## Implementation structure

- Cross-polarisation experiment in the rotating frame, followed by
- time-domain FID acquisition. The CP stage is preceded by wiping
- of the low-gamma spins and followed by FID acquisition with deco-
- upling of the high-gamma spins. Syntax:
- fid=cp_acquire_soft(spin_system,parameters,H,R,K)
- parameters.spins -working spins, a cell array of
- strings with high-gamma spin fi-
- rst, and low-gamma spin last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- channel, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `iscell()`, `cellfun()`, `isscalar()`, `any()`, `ismember()`.
