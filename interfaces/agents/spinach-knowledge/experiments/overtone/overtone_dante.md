# experiments/overtone/overtone_dante.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/overtone/overtone_dante.m) · [Spinach wiki](https://spindynamics.org/wiki/index.php?title=Overtone_dante.m)

Signature: `spectrum=overtone_dante(spin_system,parameters,H,R,K)`.

## Behaviour

The source implements an overtone DANTE pulse train followed by the frequency-domain acquisition helper `overtone_a`. It computes `ovt_frq=-2*spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi)`, forms `L=H+1i*R+1i*K`, and extends `parameters.Lx` across `spc_dim`. The source labels `Lx` as the X Zeeman operator on the quadrupolar nucleus. It sets the rotor period to `abs(1/rate)`, divides that period by `pulse_num` to get the pulse-cycle length, builds a pulse propagator and an intervening evolution propagator, combines them as `PE*PP`, and applies the combined propagator `n_periods*pulse_num` times to `rho0`.

The pulse offset in the average-Hamiltonian call is `2*pi*(ovt_frq-rf_frq)`. `rate` is documented in Hz; `pulse_dur` is in seconds and `pulse_amp` in rad/s. The code errors when `pulse_dur` exceeds the cycle length; equality is allowed by that check. Its `rate` field defines rotor-period timing, but the source has no separate sample-orientation input. This mapped function is not a REDOR implementation: the source constructs the DANTE pulse/evolution cycle and does not define a REDOR pulse pair or dephasing-difference acquisition. Any coupling evolution in `H` is supplied by the caller.

## Inputs and limits

Required fields: `spins` (one-element cell array), `spc_dim` (positive integer), `Lx`, `rf_frq` (offset in Hz), `rate` (non-zero real scalar in Hz), `npoints` (positive integer), `rho0`, `coil`, `sweep` (two-element real acquisition range in Hz), `pulse_dur` (positive duration in seconds), `pulse_amp` (rad/s), `pulse_num` (positive integer pulses per rotor period), and `n_periods` (positive integer active rotor periods). The consistency check requires `H`, `R`, and `K` to be numeric matrices with matching dimensions. The DANTE timing check additionally requires `pulse_dur <= abs(1/rate)/pulse_num`.

The returned `spectrum` is passed through from `overtone_a`; the wrapper supplies the acquisition sweep and point count but does not reshape the helper output.
