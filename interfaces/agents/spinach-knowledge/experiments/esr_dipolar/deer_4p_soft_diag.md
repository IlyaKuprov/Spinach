# experiments/esr_dipolar/deer_4p_soft_diag.m

- Signature: `deer_4p_soft_diag(spin_system,parameters)`

## Purpose

Runs four-pulse DEER diagnostics. It calls `deer_4p_soft_hole` through `powder`, apodises and zero-fills/FFT-transforms its four FIDs for pulse-diagnostic plots, then calls `deer_4p_soft_deer` through `powder` to obtain the echo stack. The wrapper plots the real, unphased stack and uses its SVD to plot echo and DEER components.

## Parameters / inputs

- parameters.pulse_frq -frequencies for the four
- pulses, Hz
- parameters.pulse_pwr -power levels for the four
- pulses, Hz
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
- parameters.offset -receiver offset for the time
- domain detection, Hz
- parameters.sweep -sweep width for time domain
- detection, Hz
- parameters.npoints -number of points in the free
- induction decay
- parameters.zerofill -length of the zero-filled FFT
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
