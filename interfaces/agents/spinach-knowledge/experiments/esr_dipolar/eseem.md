# experiments/esr_dipolar/eseem.m

- Signature: `fid=eseem(spin_system,parameters,H,R,K)`

## Purpose

ESEEM pulse sequence with ideal hard pulses.

## Parameters / inputs

- `parameters.npoints` — number of points to be computed.
- `parameters.timestep` — simulation time step, seconds.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.screen` — optional screen state; must be the Hermitian conjugate of the detection state. Defaults to `[]` if omitted.
- `parameters.pulse_op` — pulse operator.
- `H` — Hamiltonian matrix, received from context function.
- `R` — relaxation superoperator, received from context function.
- `K` — kinetics superoperator, received from context function.

## Outputs

- `fid` — time-domain signal whose Fourier transform is the ESEEM spectrum.

## Implementation structure

The function moves into the adjoint representation if needed, checks input consistency, and forms the Liouvillian as `L=H+1i*R+1i*K`. It applies a `pi/2` pulse to the initial state, evolves for the spin echo using `parameters.timestep/2` and `parameters.npoints-1`, applies a `pi` pulse, and performs refocusing evolution. Detection uses `parameters.coil` to produce `fid`.

Source: <https://spindynamics.org/wiki/index.php?title=eseem.m>