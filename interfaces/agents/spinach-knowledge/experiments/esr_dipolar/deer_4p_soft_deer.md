# experiments/esr_dipolar/deer_4p_soft_deer.m

- MATLAB implementation: [experiments/esr_dipolar/deer_4p_soft_deer.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_4p_soft_deer.m)

## Purpose

This is the four-pulse DEER/PELDOR simulation callback. It uses Fokker–Planck soft pulses and returns an echo stack as the third-pulse insertion time is varied; it is an experiment implementation, not a self-contained example. It is intended to be called with the spin system and Hamiltonian, relaxation, and kinetics matrices supplied by a Spinach context or a powder-simulation driver.

## Interface and parameters

`echo_stack=deer_4p_soft_deer(spin_system,parameters,H,R,K)`

The required `parameters` fields are:

- `pulse_frq`, `pulse_pwr`, `pulse_dur`, `pulse_phi`, and `pulse_rnk`: four values for the four soft pulses. The source documents frequencies in Hz, powers in rad/s, durations in seconds, phases in radians, and Fokker–Planck ranks.
- `offset`: a finite real scalar in Hz. The validator requires this field although the parameter list in the source header omits it; the implementation includes it in the pulse-frequency conversion.
- `p1_p2_gap` and `p2_p4_gap`: delays in seconds. The source describes these relative to the first, second, and third pulses; the third-pulse position is swept over the interval used by the code. The grumbler requires `p2_p4_gap` to be strictly greater than `p1_p2_gap` and the echo-window start `(parameters.p2_p4_gap-parameters.p1_p2_gap)+parameters.pulse_dur(3)-parameters.echo_time/2` to be non-negative.
- `p3_nsteps`: positive trajectory step count, yielding `p3_nsteps+1` third-pulse positions; `echo_time` in seconds and `echo_npts` steps, yielding `echo_npts+1` echo-window samples.
- `rho0` (initial state), `coil` (detection state), `spins` (irradiated spins, normally `{'E'}`), and `method` (soft-pulse propagation option: `'expm'`, `'expv'`, or `'evolution'`).

`H`, `R`, and `K` are same-sized matrices from the context function. The routine converts to Liouville representation when needed, and requires `sphten-liouv` or `zeeman-liouv` formalism.

## Sequence and signal

The pulse operators are formed from the raising operator for `parameters.spins{1}`. The implementation combines `H`, `R`, and `K` into the propagator generator, adjusts the four pulse frequencies for the electron reference frequency and `offset`, then applies four shaped pulses with free evolution between them. It creates a set of states at the third-pulse positions, refocuses that set, applies the fourth pulse, and samples the echo through `coil`. The trace time is referenced to the third-pulse insertion point after the second pulse.

The returned `echo_stack` has `echo_npts+1` time-sample rows and `p3_nsteps+1` third-pulse-position columns, including the initial observation and initial trajectory state. The echo window is positioned relative to the expected second-echo location. The source cautions that simulated echoes can be sharp because the simulation does not include experimental distributions, and recommends Fourier-transforming the echo before integration.

## Scope not specified by the source

The function does not prescribe pulse values, an initial-state or coil construction, a powder distribution, or an experimental receiver convention beyond its parameters. It does not define the output matrix orientation in its header or provide a stand-alone experimental example.

Source: https://spindynamics.org/wiki/index.php?title=deer_4p_soft_deer.m
