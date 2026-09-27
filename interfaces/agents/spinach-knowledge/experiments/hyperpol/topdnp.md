# experiments/hyperpol/topdnp.m

- Signature: `contact_curve=topdnp(spin_system,parameters,H,R,K)`

## Purpose

Simulate the time-optimised pulsed DNP sequence described in the cited paper and return the detected signal after each TOP DNP loop. Call this function from a powder context.

## Method

The function combines the supplied matrices as `L=H+1i*R+1i*K` and adds x-directed microwave irradiation during each pulse. Starting from `parameters.rho0`, it applies the pulse and then the delay for each loop, recording the coil observable after every loop; the initial observable is also included in the returned curve. In Hilbert space the pulse-delay propagator is precomputed and applied repeatedly; in Liouville space each step is applied in sequence. Supported formalisms are `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`.

## Parameters / inputs

- `H` — Hamiltonian matrix supplied by the context function.
- `R` — relaxation superoperator supplied by the context function.
- `K` — kinetics superoperator supplied by the context function.
- `parameters.irr_powers` — microwave amplitude (electron nutation frequency), in Hz.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.pulse_dur` — pulse duration, in seconds.
- `parameters.delay_dur` — delay duration, in seconds.
- `parameters.nloops` — number of TOP DNP loops.

## Output

- `contact_curve` — coil-state observable at the initial state and after each loop.

## Reference

- https://doi.org/10.1126/sciadv.aav6909
- Source documentation: <https://spindynamics.org/wiki/index.php?title=topdnp.m>
