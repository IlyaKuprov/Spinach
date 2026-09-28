# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1_r.m

- Signature: `xix_w_pulse_dur_ensemble_b1_r()`

## Purpose

2D parameter scan of XiX DNP in the steady state with electron-proton distance and electron Rabi frequency ensembles. Calculation time: hours.

## Physical / mathematical content

- Models an electron–proton pair at a 3.4 T W-band field and 80 K, with a trityl electron g-tensor, a proton Zeeman shift, and distance-dependent, orientation-dependent proton longitudinal relaxation. The steady-state XiX calculation detects the proton `Lz` expectation value under microwave irradiation.

## Numerical / algorithmic content

- Scans 101 microwave resonance offsets from −230 to 205 MHz and 200 electron pulse durations from 2 to 21 ns. For each of three electron–proton distances and five electron nutation frequencies, a `parfor` loop evaluates `powder(spin_system,@xixdnp_steady,localpar,'esr')` on the `rep_2ang_800pts_sph` orientation grid. The result is averaged over the B1 quadrature weights and then over the distance quadrature weights with an `r^2` Jacobian.

## Implementation structure

- Defines the two-spin system, Zeeman interactions, relaxation model, and full spin–Liouville basis separately for each distance.
- Sets the proton detection operator and XiX experiment parameters, including the inverted second-pulse phase, loop count, and shot spacing.
- Stores the offset-by-duration results across both ensembles, plots the real proton `Lz` expectation value against offset and pulse duration, and saves `xix_w_pulse_dur_ensemble_b1_r.fig`.
