# experiments/sp_acquire.m

- MATLAB source: [experiments/sp_acquire.m](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/sp_acquire.m)
- Signature: `fid=sp_acquire(spin_system,parameters,H,R,K)`

## What the routine does

This routine applies one Fokker–Planck-simulated soft pulse and then delegates time-domain signal acquisition to `acquire`. It forms `L=H+1i*R+1i*K` after conversion to the adjoint representation when needed. Pulse X and Y operators are built from `parameters.spins{1}`, extended over the Fokker–Planck spatial basis, and passed with the caller's `parameters.rho0` to `shaped_pulse_af`.

Before that call, the code replaces `parameters.pulse_frq` with `parameters.pulse_frq-parameters.offset`. The routine passes the pulse frequency, power, duration, phase, Fokker–Planck cut-off rank, and selected propagation method to `shaped_pulse_af`; it does not define an additional chirp law or a gradient-encoding sequence here.

## Required parameters and units

- `pulse_frq`: soft-pulse frequency relative to the current rotating frame, Hz.
- `pulse_phi`: pulse phase, rad; `pulse_pwr`: pulse power, rad/s; `pulse_dur`: pulse duration, s.
- `pulse_rnk`: Fokker–Planck cut-off rank; the source suggests starting at 2 and increasing until the answer stops changing.
- `offset`: transmitter/receiver offset relative to the current rotating frame, Hz.
- `sweep`: acquisition sweep width, Hz; `npoints`: positive integer number of FID samples.
- `rho0`: initial state; `coil`: detection state; `spins`: a one-element cell array of spin-label strings, whose first member defines the pulse operators.
- `method`: one of `expv`, `expm`, or `evolution` for soft-pulse propagation.
- `H`, `R`, and `K`: generators supplied by the context. Direct `sphten-liouv` and `zeeman-liouv` inputs are supported; a `zeeman-hilb` density-matrix context is also accepted because `sim2liouv` converts its generators, `rho0`, and `coil` to `zeeman-liouv` before the formalism guard.

## Output axis

`fid` is the time-domain signal returned by `acquire`. That routine uses a dwell interval of `1/sweep` seconds and `npoints-1` evolution steps, yielding `npoints` observable samples for one initial state. The function returns the signal vector, not a separately constructed time vector.

## References

- [Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/sp_acquire.m)
- [Spinach Wiki: sp_acquire.m](https://spindynamics.org/wiki/index.php?title=sp_acquire.m)
