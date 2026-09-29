# examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_field_profile_ensemble_b1_r.m)

## Purpose

Run the no-argument MATLAB function `xix_w_field_profile_ensemble_b1_r()` with Spinach and its example helpers available; it computes a steady-state XiX DNP field profile with both electron–proton distance and electron nutation-frequency (B1) averaging. The source labels its magnet section Q-band and sets `sys.magnet=3.4`; it estimates minutes of calculation time.

## Model and nested averaging

The electron–proton model uses the trityl g-tensor [2.00319 2.00319 2.00258], proton shift [0 0 5] (the source labels the proton value a ppm guess), Euler angles [0 10 0] and [0 0 10] degrees, and spin-temperature value 80. A four-point Gauss–Legendre distance quadrature uses bounds 3.5–20 Å; the B1 quadrature uses six points over 10e6–20e6 Hz. For each distance, the electron–proton coordinates are reset and `r1n_dnp` supplies distance- and orientation-dependent proton longitudinal relaxation using additional parameters `2.00230`, `1e-3`, and `52`; the source sets `r1_rates={1e3 r1n_rate}`, `r2_rates={200e3 50e3}`, `t1_t2`, `rlx_keep='diagonal'`, and `equilibrium='dibari'`. The source uses `sphten-liouv` with no approximation and a propagator chop tolerance of `1e-12`.

For every distance/B1 pair, `powder(spin_system,@xixdnp_steady,parameters,'esr')` evaluates 201 microwave offsets from −300e6 to 300e6 Hz. The XiX settings are 18 ns pulses, 10 blocks, an inverted second-pulse phase ( `pi` ), shot spacing 167e−6, and additional shift −33e6; the powder grid is `rep_2ang_800pts_sph`. The aggregation first applies normalised B1 quadrature weights, then distance weights multiplied by `r^2` and normalised by their weighted sum. The distance quadrature is explicitly labelled Å in this source; no unit is stated for the `sys.magnet` value.

## Output and dependencies

The final real DNP profile is plotted as the proton `Lz` expectation value versus microwave offset in MHz and saved to `xix_w_field_profile_ensemble_b1_r.fig`. The driver uses Spinach system/basis/detection and powder routines and the helpers `gaussleg`, `r1n_dnp`, and `xixdnp_steady`. It saves a figure, not a separate numeric data file.
