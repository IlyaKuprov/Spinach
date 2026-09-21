# experiments/holeburn.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/holeburn.m`
- Signature: `fid=holeburn(spin_system,parameters,H,R,K)`
- Total lines: 183

## Purpose

Hole burning experiment -a soft pulse follwed by a hard pi/2 observation pulse. The soft pulse is simulated using Fokker- Planck formalism. Syntax: fid=holeburn(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters.pulse_frq -frequency of the soft pulse, Hz
- parameters.pulse_phi -phase of the soft pulse, rad
- parameters.pulse_pwr -power of the soft pulse, rad/s
- parameters.pulse_dur -duration of the sof pulse, s
- parameters.pulse_rnk -Fokker-Planck cut-off rank
- parameters.offset -receiver offset for the time
- domain detection, Hz
- parameters.sweep -sweep width for time domain
- detection, Hz
- parameters.npoints -number of points in the free
- induction decay
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.method -soft puse propagation method,
- 'expv' for Krylov propagation,
- 'expm' for exponential propa-
- gation, 'evolution' for Spin-
- ach evolution function
- H -Hamiltonian matrix, received from context
- function
- R -relaxation superoperator, received from
- context function
- K -kinetics superoperator, received from
- context function
- Output:
- fid -free induction decay seen by the state speci-
- fied in parameters coil after the hole burning
- pulse followed by a hard pi/2 pulse.
- Note: the rank parameter should be increased until conver-
- gence is achieved in the output.

## Implementation structure

- Hole burning experiment -a soft pulse follwed by a hard pi/2
- observation pulse. The soft pulse is simulated using Fokker-
- Planck formalism. Syntax:
- fid=holeburn(spin_system,parameters,H,R,K)
- parameters.pulse_frq -frequency of the soft pulse, Hz
- parameters.pulse_phi -phase of the soft pulse, rad
- parameters.pulse_pwr -power of the soft pulse, rad/s
- parameters.pulse_dur -duration of the sof pulse, s
- parameters.pulse_rnk -Fokker-Planck cut-off rank
- parameters.offset -receiver offset for the time
- domain detection, Hz
- parameters.sweep -sweep width for time domain

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `speye()`, `shaped_pulse_af()`, `step()`, `acquire()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `iscell()`, `ischar()`, `isscalar()`.
