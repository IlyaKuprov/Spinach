# experiments/esr_dipolar/deer_4p_soft_deer.m

- Signature: `echo_stack=deer_4p_soft_deer(spin_system,parameters,H,R,K)`

## Purpose

Simulates the four-pulse DEER/PELDOR sequence and returns `echo_stack`, an echo sampled over the third-pulse-position interval.

## Physical / mathematical content

The function converts to Liouville representation as needed, forms `L = H + 1i*R + 1i*K`, then applies the four shaped pulses and free-evolution intervals specified by the parameters. It samples the second echo around its expected position to form the echo stack.

## Numerical / algorithmic content

Soft pulses are propagated with `shaped_pulse_af`, and free intervals with `evolution`; the propagation method is selected by `parameters.method`.

## Parameters / inputs

- parameters.pulse_frq -frequencies for the four
- pulses, Hz
- parameters.pulse_pwr -power levels for the four
- pulses, rad/s
- parameters.pulse_dur -durations for the four
- pulses, seconds
- parameters.pulse_phi -initial phases for the four
- pulses, radians
- parameters.pulse_rnk -Fokker-Planck ranks for the
- four pulses
- parameters.p1_p2_gap -time between the end of the
- first and the start of the
- second pulse, seconds
- parameters.p2_p4_gap -time between the end of the
- second the start of the third
- pulse, seconds
- parameters.p3_nsteps -number of third pulse posi-
- tions in the interval between
- the first echo and the fourth
- pulse
- parameters.echo_time -time to sample around the ex-
- pected second echo position
- parameters.echo_npts -number of points in the second
- echo discretization
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.spins -irradiated spins, normally {'E'}
- parameters.method -soft puse propagation method,
- 'expv' for Krylov propagation,
- 'expm' for exponential propa-
- gation, 'evolution' for Spin-
- ach evolution function
- H -Hamiltonian matrix, received from context function
- R -relaxation superoperator, received from context function
- K -kinetics superoperator, received from context function

## Outputs

- echo_stack -DEER echo stack, a matrix with p3_nsteps echoes
- with echo_npts points each
- Note: for the method, start with 'expm', change to 'expv' if the
- calculation runs out of memory, and use 'evolution' as the
- last resort.
- Note: simulated echoes tend to be sharp and hard to catch becau-
- se simulation does not have distributions in experimental
- parameters. Fourier transforming the echo prior to integ-
- ration is recommended.
- Note: the time in the DEER trace refers to the third pulse inser-
- tion point, after end of the second pulse.
