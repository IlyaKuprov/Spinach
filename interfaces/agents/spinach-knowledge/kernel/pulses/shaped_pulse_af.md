# kernel/pulses/shaped_pulse_af.m

- MATLAB source: [kernel/pulses/shaped_pulse_af.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/shaped_pulse_af.m)
- Spinach wiki: [shaped_pulse_af.m](https://spindynamics.org/wiki/index.php?title=shaped_pulse_af.m)
- Signature: `[rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,rf_amp_list,rf_dur_list,rf_phi,max_rank,method)`
- Reference: [Eq. 33](http://dx.doi.org/10.1016/j.jmr.2016.07.005)

## Purpose

Propagates a shaped RF pulse in amplitude-frequency coordinates using the Fokker–Planck formalism in Eq. 33 of the cited paper.

## Construction and propagation

The source sets the phase-coordinate dimension to `2*max_rank+1`, obtains phase coordinates and a derivative operator with `fourdif`, and adds `rf_phi` to the phase coordinates. It builds the background term from `L0`, the phase-dependent RF term from `cos(phases)*Lx + sin(phases)*Ly`, and the phase-turning generator from the derivative operator. For slice `n`, the generator passed to the propagator is `F0 + rf_amp_list(n)*F1 + 2i*pi*rf_frq_list(n)*M`, applied for `rf_dur_list(n)`.

The frequency, amplitude, and duration lists describe successive piecewise-constant slices and must have equal lengths. Their time discretisation must also be sufficiently fine to reproduce the intended waveform; check convergence as the slices are refined, independently of `max_rank` convergence. This implementation does not apply a separate temporal window or filter. Its supported formalisms are state-vector based: `sphten-liouv`, `zeeman-liouv`, or `zeeman-wavef`.

## Inputs and outputs

- `L0` — background drift Liouvillian; `Lx` and `Ly` — X and Y projections of the RF operator.
- `rho` — initial state vector or horizontal stack of state vectors.
- `rf_frq_list` — RF frequencies in Hz, relative to the offsets and/or rotating frames used to construct `L0`.
- `rf_amp_list` — RF amplitudes in radians per second; `rf_dur_list` — slice durations in seconds.
- `rf_phi` — phase of the first pulse slice, added to the phase coordinates.
- `max_rank` — finite positive integer truncation rank. Increase it until the result stops changing; the source describes 2 as a starting point, not a guaranteed converged setting.
- `method` — `'expv'` (Krylov propagation), `'expm'` (explicit exponential propagation), or `'evolution'` (Spinach evolution function).
- `rho` — propagated state; optional `traj` — a `1 x (nsteps+1)` cell array including the initial state; optional `P` — effective pulse propagator, available with `method='expm'` and described by the source as expensive.

The source warns that the sign convention for `rf_frq_list` must agree with the offsets and rotating frames in `L0`; it notes that the wrong sign can place the pulse far from the intended location.
