# experiments/esr_dipolar/deer_4p_soft_hole.m

- Signature: `fids=deer_4p_soft_hole(spin_system,parameters,H,R,K)`

## Purpose

Pulse diagnostics for the four-pulse DEER/PELDOR pulse sequence. This function shows how soft pulses affect the magnetisation of the sample. It is a hypothetical experiment where a user-specified soft pulse is performed, immediately followed by an ideal `pi/2` pulse on all spins and infinite-bandwidth time-domain detection.

## Physical / mathematical content

- The function composes the Liouvillian as `L=H+1i*R+1i*K`.
- Pulse operators are formed from `L+` for `parameters.spins{1}`, with `Ex=(Ep+Ep')/2` and `Ey=(Ep-Ep')/2i`.
- Pulse frequency offsets are calculated using the electron magnetisation, `parameters.pulse_frq`, and `parameters.offset`.

## Numerical / algorithmic content

- The function moves into the adjoint representation if needed and checks input consistency.
- Four soft-pulse states are calculated with `shaped_pulse_af`, one for each set of pulse frequency, power, duration, phase, and Fokker-Planck rank.
- A hard `pi/2` pulse is applied using `step(spin_system,Ey,[rho1 rho2 rho3 rho4],pi/2)`. The resulting states are passed to `acquire`.

## Parameters / inputs

- `parameters.pulse_frq` — frequencies for the four pulses, Hz.
- `parameters.pulse_pwr` — power levels for the four pulses, rad/s.
- `parameters.pulse_dur` — durations for the four pulses, seconds.
- `parameters.pulse_phi` — initial phases for the four pulses, radians.
- `parameters.pulse_rnk` — Fokker-Planck ranks for the four pulses.
- `parameters.offset` — receiver offset for time-domain detection, Hz.
- `parameters.sweep` — sweep width for time-domain detection, Hz.
- `parameters.npoints` — number of points in the free induction decay.
- `parameters.spins` — irradiated spins, normally `{'E'}`.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.method` — soft pulse propagation method: `'expv'` for Krylov propagation, `'expm'` for exponential propagation, or `'evolution'` for the Spinach evolution function.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.

## Outputs

- `fids` — four free induction decays that should be apodised and Fourier transformed.

## Implementation structure

The function applies each soft pulse to `parameters.rho0`, combines the four resulting states for the hard pulse, and acquires the resulting free induction decays. The consistency checks require `H`, `R`, and `K` to be matrices of the same dimensions and the formalism to be `sphten-liouv` or `zeeman-liouv`.

Source: <https://spindynamics.org/wiki/index.php?title=deer_4p_soft_hole.m>