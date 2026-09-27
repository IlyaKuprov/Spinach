# experiments/esr_dipolar/deer_3p_soft_diag.m

- Signature: `deer_3p_soft_diag(spin_system,parameters)`

## Purpose

Runs the three-pulse diagnostic sequence and an echo-stack calculation. It calls `deer_3p_soft_hole` through `powder`, apodises and zero-fills/FFT-transforms its four FIDs for pulse diagnostic plots, then calls `deer_3p_soft_deer` through `powder`. It plots the unphased echo stack and SVD-derived echo and DEER components against their respective axes.

## Parameters / inputs

- parameters.pulse_frq -frequencies for the three
- pulses, Hz
- parameters.pulse_pwr -power levels for the three
- pulses, Hz
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
- parameters.zerofill -length of the zero-filled FFT
- parameters.rho0 -initial state
- parameters.coil -detection state
- parameters.p1_p3_gap -time between the first and the
- third pulses, seconds
- parameters.p2_nsteps -number of second pulse posi-
- tions in the interval between
- the first and the third pulse
- parameters.echo_time -time to sample around the ex-
- pected echo position
- parameters.echo_npts -number of points in the echo
- discretization
- parameters.method -soft puse propagation method,
- 'expv' for Krylov propagation,
- 'expm' for exponential propa-
- gation, 'evolution' for Spin-
- ach evolution function
- parameters.assumptions -Hamiltonian generation assump-
- tions, use 'deer' to keep two-
- electron flip-flop terms and
- 'deer-zz' to drop them

## Outputs

- Figure 1: pulse diagnostics
- Figure 2: DEER echo stack
- Figure 3: principal components of the stack, echo
- Figure 4: principal components of the stack, DEER
- Note: for the method, start with 'expm', change to 'expv' if the
- calculation runs out of memory, and use 'evolution' as the
- last resort.
- Note: simulated echoes tend to be sharp and hard to catch becau-
- se simulation does not have distributions in experimental
- parameters. Fourier transforming the echo prior to integ-
- ration is recommended.
- Note: the time in the DEER trace refers to the second pulse inser-
- tion point, after end of first pulse.
