# examples/optimal_control/features_trapezium.m

- Signature: `features_trapezium()`

## Purpose

This example configures a piecewise-linear optimal-control pulse for transfer from proton `Lz` to fluorine `Lz` in a three-spin H–C–F system. The source describes the waveform treatment as derivatives of a Lie-group product quadrature and cites [doi:10.1016/j.jmr.2023.107478](https://doi.org/10.1016/j.jmr.2023.107478).

## Spin model and transfer

The model contains `1H`, `13C` and `19F` at 9.4 T. All chemical shifts are 0.0 ppm; the H–C and C–F scalar couplings are 140 Hz and −160 Hz. The basis is `sphten-liouv` with approximation `none`. Normalised `Lz` states on spins 1 and 3 form the initial and target states.

## Piecewise-linear control design

The six controls are `Lx` and `Ly` on each isotope, with channel map `[1;1;2;2;3;3]`. Five configured power levels run from 0.8 × 10³ × 2π to 1.2 × 10³ × 2π rad/s. The time grid has 50 slices of 0.2 ms (10 ms total), represented by 51 control values at slice endpoints. The code sets `integrator='trapezium'`, `method='lbfgs'`, `max_iter=100`, and an `SNS` penalty of weight 100. It passes a random `6 x 51` guess to `fmaxnewton(spin_system,@grape_xy,guess)` and rescales the returned waveform by the mean configured power level.

For the follow-up simulation, each slice uses generators formed from the left, midpoint and right endpoint controls, then advances the state with `step`. The script computes and reports the real initial-to-target overlap. The iteration limit and requested report do not establish convergence or provide a measured fidelity in this source.

Source: [examples/optimal_control/features_trapezium.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/features_trapezium.m).
The initial and target operator vectors are constructed with unweighted `coil_state` before norm normalisation. This single-substance model has no chemical reactions.
