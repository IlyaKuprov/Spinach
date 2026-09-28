# experiments/esr_dipolar/oopeseem.m

- Signature: `fid=oopeseem(spin_system,parameters,H,R,K)`

## Purpose

Out-of-phase ESEEM pulse sequence with the first pulse set to `pi/4` to probe two-electron correlations in the initial condition. The source comment gives the syntax as `fid=eseem(spin_system,parameters,H,R,K)`; the function declaration is `fid=oopeseem(spin_system,parameters,H,R,K)`.

## Parameters / inputs

- `parameters.npoints` — number of time points to be computed.
- `parameters.timestep` — simulation time step, seconds.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.screen` — optional screen state; must be the Hermitian conjugate of the detection state. Defaults to `[]` if absent.
- `parameters.pulse_op` — pulse operator `A`; the pulse propagators are `exp(-i*A*pi)` and `exp(-i*A*pi/4)`.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.

## Outputs

- `fid` — OOP-ESEEM time trace.

## Implementation structure

The function calls `sim2liouv()`, checks input consistency, and forms `L=H+1i*R+1i*K`. It applies a `pi/4` pulse to `parameters.rho0`, evolves the state for a spin echo, applies a `pi` pulse, and evolves it again with refocusing. Detection uses `parameters.coil` to produce `fid`.

The sequence uses ideal pulses; the source comment says to replace them with `shaped_pulse_af()` to have soft pulses instead.

[Source documentation](https://spindynamics.org/wiki/index.php?title=oopeseem.m)