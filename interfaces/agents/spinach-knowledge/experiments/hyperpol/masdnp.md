# experiments/hyperpol/masdnp.m

- Signature: `dnp=masdnp(spin_system,parameters)`

## Purpose

Simulate magic-angle-spinning dynamic nuclear polarisation and return the spherical-grid-averaged enhancement of the detected state relative to thermal equilibrium. The function uses code donated by Fred Mentink; cite his papers when using it.

## Method

For each orientation in the selected spherical grid, the function constructs the rotor-frame Hamiltonian stack and propagates through one rotor period with microwave irradiation and relaxation. It forms an effective rotor-period Liouvillian from the period propagator, evolves the thermal-equilibrium state for `parameters.mw_time`, then averages the detected magnetisation over a rotor period. The orientation-weighted enhancements are accumulated with a `parfor` loop.

## Parameters / inputs

- `parameters.spins` — spins irradiated by the microwave field.
- `parameters.rate` — spinning rate in Hz.
- `parameters.axis` — spinning-axis direction vector.
- `parameters.max_rank` — rotor-discretisation grid rank.
- `parameters.mw_pwr` — microwave power in rad/s.
- `parameters.mw_frq` — microwave frequency in Hz.
- `parameters.mw_time` — microwave irradiation duration before the average magnetisation is computed, in seconds.
- `parameters.grid` — name of the spherical-averaging grid.
- `parameters.coil` — detection state.
- `parameters.verbose` — set to 1 to enable diagnostic output.

## Output

- `dnp` — enhancement of the detected state relative to thermal equilibrium, averaged over the spherical grid and rotor period.

## Notes

- Increase the rotor rank and spherical-grid size until the answer stops changing; both may need to be large.
- Call this function directly, without a context wrapper.
- Authors: Ilya Kuprov and Fred Mentink.
- Source documentation: <https://spindynamics.org/wiki/index.php?title=masdnp.m>
