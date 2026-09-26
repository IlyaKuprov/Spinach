# experiments/sp_acquire.m

- Signature: `fid=sp_acquire(spin_system,parameters,H,R,K)`

## Purpose

Soft pulse followed by acquisition. The soft pulse is simulated using the Fokker-Planck formalism. Syntax: fid=sp_acquire(spin_system,parameters,H,R,K)

## Physical / mathematical content

- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

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
