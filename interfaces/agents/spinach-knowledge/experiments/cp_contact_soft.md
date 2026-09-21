# experiments/cp_contact_soft.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/cp_contact_soft.m`
- Signature: `contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)`
- Total lines: 150

## Purpose

Cross-polarisation experiment in the rotating frame. Applies a soft pi/2 pulse using the specified operators, then evolves the system with the specified spin-lock terms added to the Liovilli- an. The contact curve is returned. Syntax: contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.spins -working spins, a cell array of
- strings with high-gamma spins
- first and low-gamma spins last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- spins, Hz
- parameters.cp_pwr -nutation frequencies on the two
- channels during the CP contact
- time, a two-element vector, Hz
- parameters.timestep -time step of the CP contact ti-
- me, seconds
- parameters.nsteps -number of time steps to take
- during the CP contact time
- parameters.rho0 -initial state, the state of the
- low-gamma spins will be wiped
- parameters.coil -detection state vector
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function
- Output:
- contact_curve -contact curve detected on the coil
- state specified in parameters.coil

## Implementation structure

- Cross-polarisation experiment in the rotating frame. Applies a
- soft pi/2 pulse using the specified operators, then evolves the
- system with the specified spin-lock terms added to the Liovilli-
- an. The contact curve is returned. Syntax:
- contact_curve=cp_contact_soft(spin_system,parameters,H,R,K)
- parameters.spins -working spins, a cell array of
- strings with high-gamma spins
- first and low-gamma spins last,
- for example {'1H','13C'}
- parameters.hi_pwr -nutation frequency of the exci-
- tation pulse on the high-gamma
- spins, Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `decouple()`, `operator()`, `speye()`, `step()`, `evolution()`, `ismatrix()`, `all()`, `isfield()`, `iscell()`, `cellfun()`, `isscalar()`, `isrow()`, `any()`.
