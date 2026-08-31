# experiments/sp_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/sp_acquire.m`
- Signature: `fid=sp_acquire(spin_system,parameters,H,R,K)`
- Total lines: 178

## Purpose

Soft pulse followed by acquisition. The soft pulse is simulated using the Fokker-Planck formalism. Syntax: fid=sp_acquire(spin_system,parameters,H,R,K)

## Physical / mathematical content

- This file belongs to the `experiments` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 54-55: Move into adjoint representation if needed; implemented by `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 57-58: Check consistency; implemented by `grumble(spin_system,parameters,H,R,K)`.
- Lines 60-61: Compose Liouvillian; implemented by `L=H+1i*R+1i*K`.
- Lines 63-64: Pulse operators; implemented by `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 68-69: A soft pulse; implemented by `parameters.pulse_frq=parameters.pulse_frq-parameters.offset`.
- Lines 74-75: Acquisition; implemented by `fid=acquire(spin_system,parameters,H,R,K)`.

### Key state/data transformations

- Lines 55: computes `[spin_system,parameters,H,R,K]` using `[spin_system,parameters,H,R,K]=sim2liouv(spin_system,parameters,H,R,K)`.
- Lines 61: computes `L` using `L=H+1i*R+1i*K`.
- Lines 64: computes `Ep` using `Ep=operator(spin_system,'L+',parameters.spins{1})`.
- Lines 66: computes `Ex` using `Ex=(Ep+Ep')/2; Ey=(Ep-Ep')/2i`.
- Lines 69: computes `parameters.pulse_frq` using `parameters.pulse_frq=parameters.pulse_frq-parameters.offset`.
- Lines 70-73: computes `parameters.rho0` using `parameters.rho0=shaped_pulse_af(spin_system,L,Ex,Ey,parameters.rho0, parameters.pulse_frq,parameters.pulse_pwr, parameters.pulse_dur,parameters.pulse_phi, parameters.pul…`.
- Lines 75: computes `fid` using `fid=acquire(spin_system,parameters,H,R,K)`.

### Local helper functions

- Line 80: `grumble()` — `function grumble(spin_system,parameters,H,R,K)`.
  - Representative operation: `if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))|| (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.
  - Representative operation: `(~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))`.

## Parameters / inputs

- parameters.pulse_frq -frequency of the soft pulse,
- relative to the frequency of
- the current rotating frame, Hz
- parameters.pulse_phi -phase of the soft pulse, rad
- parameters.pulse_pwr -power of the soft pulse, rad/s
- parameters.pulse_dur -duration of the sof pulse, s
- parameters.pulse_rnk -Fokker-Planck cut-off rank,
- a small integer: start with 2
- and increase until the answer
- stops changing
- parameters.offset -transmitter / receiver offset
- for the time domain pulses and
- detection, relative to the cur-
- rent rotating frame, Hz
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

## Outputs

- fid -dynamics of the coil state as a function of time

## Implementation structure

- Soft pulse followed by acquisition. The soft pulse is simulated
- using the Fokker-Planck formalism. Syntax:
- fid=sp_acquire(spin_system,parameters,H,R,K)
- parameters.pulse_frq -frequency of the soft pulse,
- relative to the frequency of
- the current rotating frame, Hz
- parameters.pulse_phi -phase of the soft pulse, rad
- parameters.pulse_pwr -power of the soft pulse, rad/s
- parameters.pulse_dur -duration of the sof pulse, s
- parameters.pulse_rnk -Fokker-Planck cut-off rank,
- a small integer: start with 2
- and increase until the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sim2liouv()`, `grumble()`, `operator()`, `speye()`, `shaped_pulse_af()`, `acquire()`, `ismatrix()`, `all()`, `ismember()`, `isfield()`, `iscell()`, `ischar()`, `isscalar()`.
