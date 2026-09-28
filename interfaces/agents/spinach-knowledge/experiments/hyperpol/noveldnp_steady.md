# experiments/hyperpol/noveldnp_steady.m

- Signature: `dnp=noveldnp_steady(spin_system,parameters,H,R,K)`

## Purpose

Compute the steady-state detected observable as a function of microwave resonance offset for nuclear spin orientation via electron spin locking (NOVEL) or the pulsed solid effect (SE). Call this function from a powder context. The relaxation superoperator `R` must be thermalised to a finite temperature.

## Method

For each value in `parameters.el_offs`, the function adds the resonance offset and field-profile shift to the electron `Lz` term in the Liouvillian `H + 1i*R + 1i*K`. It constructs a pulse-period propagator for the selected protocol: NOVEL applies an x-directed flip pulse and a -y contact pulse, with an optional -x flipback pulse; SE applies the -y contact pulse only. The shot-spacing delay is then included, and `steady(...,'newton')` finds the steady state. The output is the coil-state observable for each offset.

## Parameters / inputs

- `H` — Hamiltonian matrix supplied by the context function.
- `R` — relaxation superoperator supplied by the context function; thermalise it to a finite temperature.
- `K` — kinetics superoperator supplied by the context function.
- `parameters.irr_powers` — microwave amplitude (electron nutation frequency), in Hz.
- `parameters.coil` — detection-state column vector.
- `parameters.contact_dur` — contact-pulse duration, in seconds.
- `parameters.shot_spacing` — delay between microwave-irradiation periods.
- `parameters.flippulse` — 0 for solid effect (no flip pulse), 1 for NOVEL (90-degree flip pulse).
- `parameters.flipback` — 0 for NOVEL without a flipback pulse, 1 for NOVEL with a flipback pulse.
- `parameters.pulse_dur` — flip-pulse duration, in seconds; required for NOVEL.
- `parameters.addshift` — field-profile centring shift.
- `parameters.el_offs` — array of microwave resonance offsets.

## Output

- `dnp` — steady-state observable on the detection state as a function of microwave resonance offset.

## References

- https://doi.org/10.1016/0022-2364(88)90190-4
- https://doi.org/10.1063/1.5000528
- Source documentation: <https://spindynamics.org/wiki/index.php?title=noveldnp_steady.m>
