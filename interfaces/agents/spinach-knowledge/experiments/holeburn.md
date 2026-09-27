# experiments/holeburn.m

- Signature: `fid=holeburn(spin_system,parameters,H,R,K)`

## Purpose

Simulates a hole-burning experiment: a soft pulse modeled using the Fokker–Planck formalism, followed by a hard π/2 observation pulse and acquisition of a free induction decay (FID).

## Parameters / inputs

- `parameters.pulse_frq` — soft-pulse frequency, Hz.
- `parameters.pulse_phi` — soft-pulse phase, rad.
- `parameters.pulse_pwr` — soft-pulse power, rad/s.
- `parameters.pulse_dur` — soft-pulse duration, s.
- `parameters.pulse_rnk` — Fokker–Planck cut-off rank.
- `parameters.offset` — receiver offset for time-domain detection, Hz.
- `parameters.sweep` — sweep width for time-domain detection, Hz.
- `parameters.npoints` — number of points in the FID.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.method` — soft-pulse propagation method: `'expv'` for Krylov propagation, `'expm'` for exponential propagation, or `'evolution'` for the Spinach evolution function.
- `parameters.spins` — irradiated spins, specified as a one-element cell array containing a character string.
- `parameters.spc_dim` — spatial dimension used when extending the pulse operator.
- `H` — Hamiltonian matrix received from the context function.
- `R` — relaxation superoperator received from the context function.
- `K` — kinetics superoperator received from the context function.

## Output

- `fid` — free induction decay detected by `parameters.coil` after the hole-burning soft pulse and a hard π/2 pulse.

Increase `parameters.pulse_rnk` until the output converges.

## Implementation

The function moves the inputs into the adjoint representation if needed and checks their consistency. It forms `L=H+1i*R+1i*K`, constructs pulse operators for `parameters.spins{1}`, and extends them across the spatial dimension. After subtracting `parameters.offset` from `parameters.pulse_frq`, it applies the soft pulse with `shaped_pulse_af` using the specified propagation method. It then applies a hard π/2 pulse about the y-axis and calls `acquire` to obtain the FID.

Source: <https://spindynamics.org/wiki/index.php?title=holeburn.m>

Contact: ilya.kuprov@weizmann.ac.il