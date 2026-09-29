# experiments/hyperpol/solid_effect.m

- Signature: `answer=solid_effect(spin_system,parameters)`
- Canonical MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/solid_effect.m

## Purpose and model

Implements the large-scale solid-effect DNP model cited below. The source documents a system with one electron and one nuclear spin type, potentially with many nuclei, and accepts the `sphten-liouv` or `zeeman-liouv` formalism. It constructs the `H+`, `H0`, and `H-` Hamiltonian sectors, adds the electron microwave terms and the electron/nuclear Zeeman terms, then selects an exact or average-Hamiltonian treatment. These are spin-dynamics observables; this function does not acquire an MRI image or FID. Call `solid_effect(spin_system,parameters)` directly: it constructs its own Liouvillian and is not a callback for `liquid`, `powder`, or another context wrapper.

## Inputs and numerical choices

Required parameters include `mw_pwr` (source-documented microwave power in rad/s), `nuclear_frq` (nuclear Zeeman frequency in rad/s), `theory`, and `calc_type`. The code accepts `theory='exact'`, `ah_first_order`, `ah_second_order`, `ah_third_order`, `kb_first_order`, `kb_second_order`, `kb_third_order`, and `matrix_log`; these are the exact strings checked by its grumbler. The microwave terms in the source carry a factor of 0.25. An optional `coil` supplies detection states; if omitted, the function builds longitudinal `Lz` detection channels for all spins. `time_step` (positive seconds) and positive integer `n_steps` are used for sampled modes.

The source header describes `calc_type='time_dependence'` and `calc_type='steady_state'`; the executable code also has a `calc_type='trajectory'` branch. Time dependence returns coil-observable channels through the multichannel evolution path (one row per coil, with the initial sample plus `n_steps` propagated samples). The steady-state branch returns the coil-resolved result through the evolution total-output path; the trajectory branch returns state evolution rather than a coil signal. These are calculated outputs, not measured magnetisation.

## Spatial encoding and references

There are no gradient, field-of-view, k-space, or FID parameters or outputs in this function. The source gives theory labels and units, but no example numeric parameter set or computed result; none is claimed here.

- Large-scale solid-effect formalism: https://doi.org/10.1039/C2CP23233B
- Spin Dynamics Wiki: https://spindynamics.org/wiki/index.php?title=solid_effect.m
