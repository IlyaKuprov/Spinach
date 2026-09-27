# experiments/esr_dipolar/deer_3p_soft_hole.m

- Signature: `fids=deer_3p_soft_hole(spin_system,parameters,H,R,K)`

## Purpose

Computes pulse diagnostics for the three-pulse DEER/PELDOR sequence. It evaluates each specified soft pulse from `parameters.rho0`, applies a common ideal `pi/2` hard pulse about `Ey` to the reference and three responses, then acquires the four FIDs.

## Physical / mathematical content

The function converts to Liouville representation as needed and forms `L = H + 1i*R + 1i*K`. Each response is generated with `shaped_pulse_af` using its pulse frequency, power, duration, phase, Fokker–Planck rank, and selected method, followed by time-domain acquisition with `acquire`.

## Numerical / algorithmic content

The returned FIDs are intended for pulse diagnostics; the source recommends apodising and Fourier transforming them. This routine does not calculate the DEER echo stack.

## Parameters / inputs

- parameters.pulse_frq -frequencies for the three
- pulses, Hz
- parameters.pulse_pwr -power levels for the three
- pulses, rad/s
- parameters.pulse_dur -durations for the three
- pulses, seconds
- parameters.pulse_phi -initial phases for the three
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the
- three pulses
- parameters.offset -receiver offset for the time
- domain detection, Hz
- parameters.sweep -sweep width for time domain
- detection, Hz
- parameters.npoints -number of points in the free
- induction decay
- parameters.spins -irradiated spins, normally {'E'}
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.method -soft puse propagation method,
- 'expv' for Krylov propagation,
- 'expm' for exponential propa-
- gation, 'evolution' for Spin-
- ach evolution function
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- fids -three free induction decays that should be apo-
- dised and Fourier transformed
- Note: for the method, start with 'expm', change to 'expv' if the
- calculation runs out of memory, and use 'evolution' as the
- last resort.
