# experiments/hyperpol/beamdnp.m

- Signature: `contact_curve=beamdnp(spin_system,parameters,H,R,K)`

## Purpose

Simulates the BEAM DNP experiment described in [Science Advances](https://doi.org/10.1126/sciadv.abq0536). Called from a powder context. See also the [Spinach documentation](https://spindynamics.org/wiki/index.php?title=beamdnp.m).

## Physical / mathematical content

The routine applies an electron 90-degree flip pulse under microwave irradiation along Y, then repeats a two-pulse BEAM DNP block with alternating microwave irradiation along +X and −X. It records the detection-state overlap before the blocks and after each block.

## Numerical / algorithmic content

The Liouvillian is `L=H+1i*R+1i*K`. Electron control operators are constructed from the electron raising operator; the flip-pulse duration is `1/(4*parameters.irr_powers)`. The two block Liouvillians add and subtract `2*pi*parameters.irr_powers*Ex`, respectively. In Hilbert space, a complete-block propagator is precomputed and applied as `P*rho*P'`; in Liouville space, the two event propagators are precomputed and applied in sequence. Other formalisms produce an error.

## Parameters / inputs

- `spin_system` — spin system supplied by the calling context.
- `H` — Hamiltonian matrix, received from the context function.
- `R` — relaxation superoperator, received from the context function.
- `K` — kinetics superoperator, received from the context function.
- `parameters.irr_powers` — microwave amplitude (electron nutation frequency), Hz; a positive real scalar.
- `parameters.rho0` — initial state.
- `parameters.coil` — detection state.
- `parameters.pulse_dur` — two pulse durations in seconds, supplied as a row vector of positive elements.
- `parameters.nloops` — number of BEAM DNP blocks; a positive integer.

## Output

- `contact_curve` — time dependence of the coil-state overlap, containing the initial value and one value after each BEAM DNP block (`parameters.nloops+1` values).

## Implementation structure

The routine checks its inputs, builds the electron control operators, applies the initial 90-degree pulse, and then runs the requested blocks while recording `hdot(parameters.coil,rho)`. It supports `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv` formalisms.

Source contacts: venkata-subbarao.redrouthu@uni-konstanz.de; ilya.kuprov@weizmann.ac.il; guinevere.mathies@uni-konstanz.de.