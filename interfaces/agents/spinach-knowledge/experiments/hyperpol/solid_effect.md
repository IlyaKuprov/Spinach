# experiments/hyperpol/solid_effect.m

- Signature: `answer=solid_effect(spin_system,parameters)`

## Purpose

Simulate solid-effect dynamic nuclear polarisation with the large-scale formalism described in the cited paper. The model requires exactly one electron and one nuclear spin type, but can include many nuclei.

## Method

The function constructs the `H+`, `H0`, and `H-` Hamiltonian components, includes the microwave terms, and combines them using either the exact electron rotating-frame Hamiltonian or the selected average-Hamiltonian theory. It starts from thermal equilibrium and computes the requested output: a time-dependent coil signal, a steady-state signal, or a trajectory. It builds its own Liouvillian and must be called directly, without a context wrapper.

## Parameters / inputs

- `parameters.mw_pwr` — microwave power, in rad/s.
- `parameters.theory` — `exact`, `ah_first_order`, `ah_second_order`, `ah_third_order`, `kb_first_order`, `kb_second_order`, `kb_third_order`, or `matrix_log`. See `average.m` for the average-Hamiltonian options.
- `parameters.nuclear_frq` — nuclear Zeeman frequency, in rad/s.
- `parameters.calc_type` — `time_dependence`, `steady_state`, or `trajectory`.
- `parameters.time_step` — time step in seconds; required for `time_dependence` and `trajectory`.
- `parameters.n_steps` — number of steps; required for `time_dependence` and `trajectory`.
- `parameters.coil` — optional detection state(s); when omitted, the function uses each spin's longitudinal `Lz` state.

## Output

- `answer` — for `time_dependence`, observables detected by the coil states at each time point; for `steady_state`, the asymptotic detected values; for `trajectory`, the evolved state trajectory.

## References

- Large-scale solid-effect formalism: https://doi.org/10.1039/C2CP23233B
- Source documentation: <https://spindynamics.org/wiki/index.php?title=solid_effect.m>
