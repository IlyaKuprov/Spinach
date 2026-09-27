# experiments/esr_dipolar/deer_3p_soft_deer.m

- Signature: `echo_stack=deer_3p_soft_deer(spin_system,parameters,H,R,K)`

## Purpose and sequence

Simulate a three-pulse DEER/PELDOR sequence with soft pulses using the Fokker-Planck formalism. The function propagates the initial state through the first pulse, samples the second-pulse position across the first-to-third-pulse gap, applies the third pulse, and records an echo window for each sampled position.

## Parameters and output

- `parameters.pulse_frq`, `parameters.pulse_pwr`, `parameters.pulse_dur`, `parameters.pulse_phi` and `parameters.pulse_rnk` give the three pulses' frequencies, powers, durations, phases and Fokker-Planck ranks. Each is a three-element vector; powers and durations must be positive.
- `parameters.p1_p3_gap` sets the first-to-third-pulse gap; `parameters.p2_nsteps` sets the number of sampled second-pulse positions.
- `parameters.echo_time` and `parameters.echo_npts` specify the echo sampling window and its number of points.
- `parameters.rho0` and `parameters.coil` are the initial and detection states; `parameters.spins` selects the irradiated spin, `parameters.offset` gives the receiver offset, and `parameters.method` selects the propagation method (`expv`, `expm` or `evolution`).
- Returns `echo_stack`, with one `parameters.echo_npts`-sample trace for each of the `parameters.p2_nsteps` positions.

## Requirements

The function is available only in Liouville space. `H`, `R` and `K` must be same-sized matrices. `parameters.spins` is a one-element cell array containing a character string; `parameters.p1_p3_gap` must be positive and `parameters.p2_nsteps` and `parameters.echo_npts` positive integers.

## Notes

- The DEER-trace time refers to the second-pulse insertion point, after the first pulse ends.
- Simulated echoes can be sharp because the simulation lacks the experimental parameter distributions; Fourier-transform the echo before integration.
- For propagation, start with `expm`, switch to `expv` if memory runs out, and use `evolution` only as a last resort.
