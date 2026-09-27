# experiments/hyperpol/noveldnp.m

- Signature: `contact_curve=noveldnp(spin_system,parameters,H,R,K)`

## Purpose

Calculate the contact-time signal for nuclear spin orientation via electron spin locking (NOVEL) or the pulsed solid effect. The `parameters.flippulse` switch selects the protocol: 0 uses the initial state directly (solid effect); 1 applies a 90-degree electron flip pulse before the contact period (NOVEL).

## Method

The function forms the Liouvillian `H + 1i*R + 1i*K` and constructs the electron transverse operators. For NOVEL it evolves the initial state under an x-directed microwave pulse for `parameters.pulse_dur`; otherwise it starts from `parameters.rho0`. It then applies a continuous spin lock along -y at the specified microwave amplitude and records the coil observable over the requested contact-time steps.

## Parameters / inputs

- `H` — Hamiltonian matrix supplied by the context function.
- `R` — relaxation superoperator supplied by the context function.
- `K` — kinetics superoperator supplied by the context function.
- `parameters.irr_powers` — microwave amplitude (electron nutation frequency), in Hz.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.timestep` — contact-curve time step, in seconds.
- `parameters.nsteps` — number of time steps in the contact curve.
- `parameters.flippulse` — 0 for solid effect (no flip pulse), 1 for NOVEL (90-degree flip pulse).
- `parameters.pulse_dur` — flip-pulse duration in seconds; required when `parameters.flippulse` is 1.

## Output

- `contact_curve` — time dependence of the coil-state observable.

## References

- https://doi.org/10.1016/0022-2364(88)90190-4
- https://doi.org/10.1063/1.5000528
- Source documentation: <https://spindynamics.org/wiki/index.php?title=noveldnp.m>
