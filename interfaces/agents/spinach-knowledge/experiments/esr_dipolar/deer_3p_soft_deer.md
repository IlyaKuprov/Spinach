# experiments/esr_dipolar/deer_3p_soft_deer.m

This function simulates a three-pulse DEER/PELDOR echo stack with soft pulses propagated by the Fokker-Planck formalism. It accepts context-supplied Hamiltonian, relaxation, and kinetics matrices, converts to Liouville representation when needed, and forms `L=H+1i*R+1i*K`. The source admits the `sphten-liouv` and `zeeman-liouv` formalisms.

## Timing and detection

After the first shaped pulse, the routine propagates a trajectory across `parameters.p1_p3_gap` at `parameters.p2_nsteps` insertion positions, applies the second shaped pulse at those positions, refocuses the trajectory, and applies the third shaped pulse. It then propagates to the echo-window start and samples `parameters.coil` over `parameters.echo_time` with `parameters.echo_npts` points. The time coordinate denotes the second-pulse insertion point after the first pulse ends. This source defines three-pulse DEER, not four-pulse DEER or CPMG/CP.

Required fields are `parameters.pulse_frq`, `parameters.pulse_pwr`, `parameters.pulse_dur`, `parameters.pulse_phi`, and `parameters.pulse_rnk` (three values each), plus `parameters.p1_p3_gap`, `parameters.p2_nsteps`, `parameters.echo_time`, `parameters.echo_npts`, `parameters.rho0`, `parameters.coil`, `parameters.spins`, `parameters.offset`, and `parameters.method`. The context arguments `H`, `R`, and `K` must be same-sized matrices. Units are Hz for pulse frequencies and receiver offset, rad/s for pulse power, seconds for durations and the gap/window, and radians for phases; ranks are integer Fokker-Planck ranks. The spin cell normally contains E; the pulse operators are built for its first entry. Methods are `expm`, `expv`, or `evolution`.

`echo_stack` contains `parameters.p2_nsteps+1` echoes with `parameters.echo_npts+1` samples per echo: trajectory propagation includes the initial state, and observable sampling includes the initial observation. The source suggests `expm`, then `expv` if memory is exhausted, then `evolution`; this is source guidance, not a performance guarantee. It also notes that simulated echoes may be narrow without experimental-parameter distributions and recommends Fourier-transforming the echo before integration.

Source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_soft_deer.m
