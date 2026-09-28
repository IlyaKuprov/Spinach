# kernel/pulses/shaped_pulse_af.m

- Signature: `[rho,traj,P]=shaped_pulse_af(spin_system,L0,Lx,Ly,rho,rf_frq_list,rf_amp_list,rf_dur_list,rf_phi,max_rank,method)`

## Purpose

Propagates a shaped RF pulse in amplitude-frequency coordinates using the Fokker–Planck formalism (Eq. 33 in the cited paper).

## Algorithm

The pulse is treated as piecewise-constant over the supplied slices. The function builds the Fokker–Planck frequency grid up to `max_rank`, combines the background Liouvillian with the RF operators and the frequency term for each slice, and propagates the state for that slice's duration. The default `'expv'` method uses Krylov propagation; `'expm'` forms explicit propagators and can return the effective pulse propagator; `'evolution'` calls Spinach's `evolution` function. The frequency list is used relative to the offsets or rotating frames included in `L0`; use the correct frequency sign, because a wrong sign can place the pulse far from its intended location.

## Parameters / inputs

- `spin_system` — Spinach system description; the source accepts `sphten-liouv`, `zeeman-liouv`, and `zeeman-wavef` formalisms.
- Input checks require `L0`, `Lx`, and `Ly` to be square and the same size, with `L0` compatible with `rho`; the frequency, amplitude, and duration lists must have equal element counts and contain finite real values, with durations strictly positive. `rf_phi` must be a finite real scalar, and `max_rank` a finite positive integer.
- `L0` — background drift Liouvillian.
- `Lx`, `Ly` — X and Y projections of the RF operator.
- `rho` — initial state vector, or a horizontal stack of state vectors.
- `rf_frq_list` — RF frequencies for the time slices, in Hz, relative to the offsets and/or rotating frames used to construct `L0`.
- `rf_amp_list` — RF amplitudes for the slices, in radians per second.
- `rf_dur_list` — slice durations, in seconds.
- `rf_phi` — phase offset applied to the pulse phases.
- `max_rank` — maximum Fokker–Planck rank; increase it until the answer stops changing (the source suggests 2 as a starting point).
- `method` — `'expv'` (default, Krylov), `'expm'` (explicit exponential), or `'evolution'` (Spinach evolution function).

## Outputs

- `rho` — final state vector or stack of state vectors.
- `traj` — when requested, a `1 x (nsteps+1)` cell array of states, including the initial state as its first entry.
- `P` — effective pulse propagator, available only with `method='expm'`; the source notes that it is expensive.

The pulse is assumed piecewise-constant, so choose a sufficiently fine time discretisation to reproduce the waveform.

For the amplitude-frequency treatment, see [Eq. 33 in the cited paper](http://dx.doi.org/10.1016/j.jmr.2016.07.005) and the [Spinach wiki page](https://spindynamics.org/wiki/index.php?title=shaped_pulse_af.m).
