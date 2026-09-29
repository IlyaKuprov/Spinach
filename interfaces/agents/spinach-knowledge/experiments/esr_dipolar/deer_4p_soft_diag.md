# experiments/esr_dipolar/deer_4p_soft_diag.m

- MATLAB implementation: [experiments/esr_dipolar/deer_4p_soft_diag.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_4p_soft_diag.m)

## Purpose

This is a four-pulse DEER/PELDOR diagnostic driver rather than a signal-returning sequence function. It first examines the frequency-selective effects of the four soft pulses, then calculates and plots a powder-averaged echo stack and its principal components.

## Interface and parameters

`deer_4p_soft_diag(spin_system,parameters)` has no declared output. Its parameter structure supplies the four-pulse settings (`pulse_frq`, `pulse_pwr`, `pulse_dur`, `pulse_phi`, `pulse_rnk`), delays (`p1_p2_gap`, `p2_p4_gap`), insertion count (`p3_nsteps`), echo window (`echo_time`, `echo_npts`), initial state (`rho0`), detection state (`coil`), one irradiated-spin selector (`spins`, a one-element cell such as `{'E'}`), receiver offset (`offset`), time-domain sweep width (`sweep`), FID point count (`npoints`), zero-filled FFT length (`zerofill`), shaped-pulse `method` (`'expm'`, `'expv'`, or `'evolution'`), and Hamiltonian-generation `assumptions`. The source gives frequencies and offsets/sweep in Hz, pulse durations in seconds, phases in radians, and echo time in seconds. It documents `echo_npts` and `p3_nsteps` as at least two.

Despite the wrapper header’s Hz label, `pulse_pwr` is consumed in rad/s: the driver passes it unchanged through both callbacks to `shaped_pulse_af`. The allowed `assumptions` labels are `'deer'` (retain two-electron flip-flop terms) and `'deer-zz'` (drop them).

## Diagnostic interpretation

The first calculation uses `deer_4p_soft_hole` under `powder`, apodises the four resulting FIDs with the `crisp` window, and Fourier-transforms them to compare the four pulse responses. The second uses `deer_4p_soft_deer` under `powder` to obtain the echo stack. The driver plots the unphased stack and uses singular-value decomposition to form its echo- and DEER-domain principal-component plots. Its documented outputs are four figures: pulse diagnostics, echo stack, echo principal components, and DEER principal components.

## Scope not specified by the source

This function does not return the FIDs, echo stack, or principal-component arrays to its caller; its declared interface is plotting-only. The source supplies no default parameter values, sample data, or numeric results, and does not define the matrix orientation of the internal stack in the diagnostic text.

Source: https://spindynamics.org/wiki/index.php?title=deer_4p_soft_diag.m
