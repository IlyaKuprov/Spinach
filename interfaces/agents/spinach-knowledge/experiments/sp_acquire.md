# experiments/sp_acquire.m

- Signature: `fid=sp_acquire(spin_system,parameters,H,R,K)`

## Purpose

Applies a soft pulse and then acquires the free induction decay. The soft pulse is simulated using the Fokker–Planck formalism.

## Numerical / algorithmic content

The routine moves the system into the adjoint representation if needed, forms `L=H+1i*R+1i*K`, and constructs the X and Y pulse operators for `parameters.spins{1}`. It subtracts `parameters.offset` from `parameters.pulse_frq`, passes the adjusted frequency and pulse settings to `shaped_pulse_af`, then calls `acquire` to produce the FID.

## Parameters / inputs

- `parameters.pulse_frq` — soft-pulse frequency relative to the current rotating frame, Hz
- `parameters.pulse_phi` — soft-pulse phase, rad
- `parameters.pulse_pwr` — soft-pulse power, rad/s
- `parameters.pulse_dur` — soft-pulse duration, s
- `parameters.pulse_rnk` — Fokker–Planck cut-off rank; start with 2 and increase until the answer stops changing
- `parameters.offset` — transmitter/receiver offset
- `parameters.sweep` — sweep width for time-domain detection, Hz
- `parameters.npoints` — number of points in the FID
- `parameters.rho0` — initial state
- `parameters.coil` — detection state
- `parameters.method` — soft-pulse propagation method: `expv`, `expm`, or `evolution`
- `parameters.spins` — working spins; the pulse operators use the first specified spin
- `H`, `R`, `K` — Hamiltonian, relaxation, and kinetics matrices, respectively, received from the context function

## Outputs

- `fid` — dynamics of the coil state as a function of time

## Reference

- [Spinach documentation](https://spindynamics.org/wiki/index.php?title=sp_acquire.m)
